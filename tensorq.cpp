// coding for bulk / tensor q for high-order

// reset viscous forces and shock heating

  for(int idim=0;idim<ndims;idim++){Fv.at(idim)=vector<double> (Fv.at(idim).size(),0.0);}
  eshock=vector<double>(eshock.size(),0.0);

// assemble force matrix to connect thermodynamic/kinematic spaces, this can be used as rhs of both energy and momentum equations

  timers.Start(TIMER_FORCE);

  for(long k=0;k<nzeroes;k++){F.at(k)=0.0;}

  for(int i=0,k=0;i<M.NCells();i++){

// update jacobian

    jacobian(i,xinit,M,S,detJ0,detDJ0);
    jacobian(i,x1,M,S,detJ,detDJ);

// evaluate energy at each integration point

    for(int gi=0;gi<S.ngi();gi++){
      double egi(0.0);
      for(int jloc=0;jloc<T.nloc();jloc++){
        egi+=e1.at(M.GlobalNode_DFEM(i,jloc))*T.value(jloc,gi);
      }
      egi=max(ECUT,egi);
    }

// update quadrature data

    for(int gi=0;gi<S.ngi();gi++){
      l.at(gi)=sqrt(V1.at(i))/S.order();
      d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);
      p.at(gi)=P(d.at(gi),egi,gamma.at(mat.at(i)-1));
      c.at(gi)=M.UpdateSoundSpeed(gamma.at(mat.at(i)-1),p.at(gi),d.at(gi));
    }

// construct force terms

    for(int idim=0;idim<M.NDims();idim++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<T.nloc();jloc++,k++){
          for(int gi=0;gi<S.ngi();gi++){
            F.at(k)+=p.at(gi)*detDJ[idim][iloc][gi]*T.value(jloc,gi)*detJ[gi]*S.wgt(gi);
          }
        }
      }
    }

  }

  timers.Stop(TIMER_FORCE);

// artificial viscosity term

  timers.Start(TIMER_VISCOSITY);

  q=vector<double> (S.ngi(),0.0);

  if(tensorq){

// tensor q

    for(int i=0;i<M.NCells();i++){

// update jacobian

      jacobian(i,x1,M,S,detJ,detDJ);

// unscaled stiffness matrix - a finite element discretisation of the grad-dot-grad operator

      vector<vector<double> > Sz(S.nloc(),vector<double> (S.nloc(),0.0));

      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<S.nloc();jloc++){
          double sij(0.0);
          for(int gi=0;gi<S.ngi();gi++){
            sij+=(detDJ[0][iloc][gi]*detDJ[0][jloc][gi]+detDJ[1][iloc][gi]*detDJ[1][jloc][gi])*detJ.at(gi)*S.wgt(gi);
          }
          Sz.at(iloc).at(jloc)=sij;
        }
      }

// viscous and corner forces (viscous force is aggregated, corner force is not)

      for(int idim=0;idim<ndims;idim++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          double f0zi(0.0);
          for(int jloc=0;jloc<S.nloc();jloc++){
            f0zi+=Sz.at(iloc).at(jloc)*u1.at(idim).at(M.GlobalNode_CFEM(i,jloc));
          }
          Fv.at(idim).at(M.GlobalNode_CFEM(i,iloc))+=f0zi;
          Fc.at(idim).at(i).at(iloc)=f0zi;
        }
      }

    }

// we need a second loop over cell to ensure we have accumulated viscous forces Fv to the global node space
// the second loop constructs the smoothness sensor on a reduction of Fc
// we also need to construct the compression/vorticity switches
// Note: we need a temporary viscous force Fvtmp as Fv is already populated and we don't want to accumulate
// on to what is in there already from the above block of code

    vector<vector<double> > Fvtmp(M.NDims(),vector<double> (nknodes,0.0));

    for(int i=0;i<M.NCells();i++){

// update jacobian

      jacobian(i,xinit,M,S,detJ0,detDJ0);
      jacobian(i,x1,M,S,detJ,detDJ);

// evaluate energy at each integration point

      for(int gi=0;gi<S.ngi();gi++){
        double egi(0.0);
        for(int jloc=0;jloc<T.nloc();jloc++){
          egi+=e1.at(M.GlobalNode_DFEM(i,jloc))*T.value(jloc,gi);
        }
        egi=max(ECUT,egi);
      }

// update quadrature data

      for(int gi=0;gi<S.ngi();gi++){
        l.at(gi)=sqrt(V1.at(i))/S.order();
        d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);
        p.at(gi)=P(d.at(gi),egi,gamma.at(mat.at(i)-1));
        c.at(gi)=M.UpdateSoundSpeed(gamma.at(mat.at(i)-1),p.at(gi),d.at(gi));
      }

// velocity at each quadrature point

      vector<vector<double> > ugi(M.NDims(),vector<double> (S.ngi(),0.0));
      for(int idim=0;idim<ndims;idim++){
        for(int gi=0;gi<S.ngi();gi++){
          for(int iloc=0;iloc<S.nloc();iloc++){
            ugi.at(idim).at(gi)+=u1.at(idim).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
          }
        }
      }

// calculate divergence field at the quadrature point

      vector<double> Cz(S.ngi(),0.0);
      for(int gi=0;gi<S.ngi();gi++){
        for(int idim=0;idim<ndims;idim++){
          for(int iloc=0;iloc<S.nloc();iloc++){
            Cz.at(gi)+=ugi.at(idim).at(gi)*detDJ.at(idim).at(iloc).at(gi);
          }
        }
      }

// curl of the velocity field at the quadrature point

      vector<double> Vz(S.ngi(),0.0);
      for(int gi=0;gi<S.ngi();gi++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          Vz.at(gi)+=detDJ.at(0).at(iloc).at(gi)*ugi.at(0).at(gi);
          Vz.at(gi)-=detDJ.at(1).at(iloc).at(gi)*ugi.at(1).at(gi);
        }
        Vz.at(gi)=abs(Vz.at(gi));
      }

// smoothness sensor

      vector<double> gp(S.ngi());
      vector<vector<double> > Fvgi(M.NDims(),vector<double> (S.ngi(),0.0));
      for(int gi=0;gi<S.ngi();gi++){
        gp.at(gi)=(V1.at(i)*c.at(gi)/(l.at(gi)*l.at(gi)));
        for(int idim=0;idim<ndims;idim++){
          for(int iloc=0;iloc<S.nloc();iloc++){
            Fvgi.at(idim).at(gi)+=Fv.at(idim).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
          }
        }
      }

      vector<double> psi0(S.ngi());
      double alpha0(0.005);
      for(int gi=0;gi<S.ngi();gi++){
        double fpnorm(sqrt(Fvgi.at(0).at(gi)*Fvgi.at(0).at(gi)+Fvgi.at(1).at(gi)*Fvgi.at(1).at(gi)));
        psi0.at(gi)=(1.0-exp(-fpnorm/(alpha0*abs(gp.at(gi)))));
      }
      double psi0max(*max_element(psi0.begin(),psi0.end()));

// compression switch

      vector<double> psi1(S.ngi());
      for(int gi=0;gi<S.ngi();gi++){
        psi1.at(gi)=((Cz.at(gi)<0.0)?1.0:0.0);
      }

// vorticity switch

      vector<double> psi2(S.ngi());
      for(int gi=0;gi<S.ngi();gi++){
        double alpha2(1.0),tmp(Vz.at(gi)/max(abs(Cz.at(gi)),1.0e-10));
        double psi2.at(gi)=(1.0/(1.0+alpha2*tmp));
      }

// viscous coefficent

      for(int gi=0;gi<S.ngi();gi++){
        mu.at(gi)=(psi0max*psi1.at(gi)*d1.at(gi)*l.at(gi)*(cq*l.at(gi)*abs(Cz.at(gi))+psi2.at(gi)*cl*c.at(gi)));
        qz.at(gi)=(mu.at(gi)*abs(Cz.at(gi))); // scalar coefficient, see eqn (42)
      }

// scaled stiffness matrix

      vector<vector<double> > Sz(S.nloc(),vector<double> (S.nloc(),0.0));

      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<S.nloc();jloc++){
          double sij(0.0);
          for(int gi=0;gi<S.ngi();gi++){
            sij+=mu.at(gi)*(detDJ[0][iloc][gi]*detDJ[0][jloc][gi]+detDJ[1][iloc][gi]*detDJ[1][jloc][gi])*detJ.at(gi)*S.wgt(gi);
          }
          Sz.at(iloc).at(jloc)=sij;
        }
      }

// viscous forces

      for(int idim=0;idim<ndims;idim++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          double f0zi(0.0);
          for(int jloc=0;jloc<S.nloc();jloc++){
            f0zi+=Sz.at(iloc).at(jloc)*u1.at(idim).at(M.GlobalNode_CFEM(i,jloc));
          }
          Fvtmp.at(idim).at(M.GlobalNode_CFEM(i,iloc))+=f0zi;
        }
      }

// shock heating

      vector<double> ek(S.nloc(),0.0),et(T.nloc(),0.0);
      for(int idim=0;idim<ndims;idim++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          ek.at(iloc)=Fc.at(idim).at(i).at(iloc)*u1.at(idim).at(M.GlobalNode_CFEM(i,iloc))*dt;
        }
      }

      S.prolongate(ek,et,T.order());

      for(int iloc=0;iloc<T.nloc();iloc++){
        eshock.at(M.GlobalNode_DFEM(i,iloc))=et.at(iloc);
      }

    }

// store viscous forces

    for(int idim=0;idim<ndims;idim++){
      for(long i=0;i<nknodes;i++){
        Fv.at(idim).at(i)=Fvtmp.at(idim).at(i);
      }
    }

  }else{

// bulk q

    for(int i=0,k=0;i<M.NCells();i++){

// update jacobian

      jacobian(i,xinit,M,S,detJ0,detDJ0);
      jacobian(i,x1,M,S,detJ,detDJ);

// evaluate energy at each integration point

      for(int gi=0;gi<S.ngi();gi++){
        double egi(0.0);
        for(int jloc=0;jloc<T.nloc();jloc++){
          egi+=e1.at(M.GlobalNode_DFEM(i,jloc))*T.value(jloc,gi);
        }
        egi=max(ECUT,egi);
      }

// update quadrature data

      for(int gi=0;gi<S.ngi();gi++){
        l.at(gi)=sqrt(V1.at(i))/S.order();
        d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);
        p.at(gi)=P(d.at(gi),egi,gamma.at(mat.at(i)-1));
        c.at(gi)=M.UpdateSoundSpeed(gamma.at(mat.at(i)-1),p.at(gi),d.at(gi));
      }

// viscous coefficient

      for(int gi=0;gi<S.ngi();gi++){
        double divu(0.0);
        for(int idim=0;idim<M.NDims();idim++){
          for(int jloc=0;jloc<S.nloc();jloc++){
            divu+=detDJ[idim][jloc][gi]*u1.at(idim).at(M.GlobalNode_CFEM(i,jloc));
          }
        }
        q.at(gi)=M.UpdateQ(l.at(gi),d.at(gi),c.at(gi),cq,cl,divu);
      }

// add on viscous forces

      for(int idim=0;idim<M.NDims();idim++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          for(int jloc=0;jloc<T.nloc();jloc++,k++){
            for(int gi=0;gi<S.ngi();gi++){
              F.at(k)+=q.at(gi)*detDJ[idim][iloc][gi]*T.value(jloc,gi)*detJ[gi]*S.wgt(gi);
            }
          }
        }
      }

    }

  }

  timers.Stop(TIMER_VISCOSITY);

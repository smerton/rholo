// Function definitions for members of the input class

// Author S. R. Merton

#include <iostream>
#include <string>
#include "input.h"
#include <fstream>

using namespace std;

// Constructor to instantiate a new input from data in the input file

Input::Input(char* inputfile){

  cout<<"Input::Input(): Loading an input from file "<<inputfile<<endl;

// read the file

  mKeyword.empty();
  mzNoh=false;
  mzSedov=false;
  mzSaltzmann=false;

  ifstream file(inputfile);

  while(mKeyword.compare("End")!=0){
    file>>mKeyword;
    cout<<"Input::Input(): input string= "<<mKeyword<<endl;
    if(mKeyword.compare("Title")==0){getline(file,mTitle);}
    if(mKeyword.compare("Description")==0){getline(file,mbuffer);mDescription.push_back(mbuffer);}
    if(mKeyword.compare("Nodes")==0){file>>mNodes[0]>>mNodes[1]>>mNodes[2];mNCells=(mNodes[0]-1)*(mNodes[1]-1);}
    if(mKeyword.compare("Domain")==0){file>>mKeyword>>mXMin>>mXMax>>mKeyword>>mYMin>>mYMax>>mKeyword>>mZMin>>mZMax;}
    if(mKeyword.compare("Boundaries")==0){for(int i=0;i<6;i++){file>>mBoundaries[i];}}
    if(mKeyword.compare("Pk")==0){file>>mPk;}
    if(mKeyword.compare("Pt")==0){file>>mPt;}
    if(mKeyword.compare("Ambient")==0){file>>mAmbient;}
    if(mKeyword.compare("Noh")==0){mzNoh=true;}
    if(mKeyword.compare("Sedov")==0){mzSedov=true;}
    if(mKeyword.compare("Saltzmann")==0){mzSaltzmann=true;}
    if(mKeyword.compare("Material")==0){
      int matno;double d,p,u[3],x[2],y[2],z[2];
      file>>matno;mMaterial.push_back(matno);
      file>>mKeyword;if(mKeyword.compare("density")==0){file>>d;}
      file>>mKeyword;if(mKeyword.compare("pressure")==0){file>>p;}
      file>>mKeyword;if(mKeyword.compare("velocity")==0){for(int i=0;i<3;i++){file>>u[i];}}
      file>>mKeyword;if(mKeyword.compare("x")==0){for(int i=0;i<2;i++){file>>x[i];}}
      file>>mKeyword;if(mKeyword.compare("y")==0){for(int i=0;i<2;i++){file>>y[i];}}
      file>>mKeyword;if(mKeyword.compare("z")==0){for(int i=0;i<2;i++){file>>z[i];}}
      mDensity.push_back(d);
      mPressure.push_back(p);
      for(int i=0;i<3;i++){mVelocity[i].push_back(u[i]);}
      mRange[0].push_back(x[0]);mRange[1].push_back(x[1]);
      mRange[2].push_back(y[0]);mRange[3].push_back(y[1]);
      mRange[4].push_back(z[0]);mRange[5].push_back(z[1]);
    }

  }

// echo the input

  cout<<endl;
  cout<<"Input::Input(): Interpreted Input:"<<endl;
  cout<<"Input::Input(): title "<<Title()<<endl;
  for(int i=0;i<NDescriptions();i++){
    cout<<"Input::Input(): description "<<Description(i)<<endl;
  }
  cout<<"Input::Input(): nodes (x-axis) "<<Nodes(0)<<endl;
  cout<<"Input::Input(): nodes (y-axis) "<<Nodes(1)<<endl;
  cout<<"Input::Input(): nodes (z-axis) "<<Nodes(2)<<endl;
  cout<<"Input::Input(): mesh (xmin,xmax) "<<XMin()<<","<<XMax()<<endl;
  cout<<"Input::Input(): mesh (ymin,ymax) "<<YMin()<<","<<YMax()<<endl;
  cout<<"Input::Input(): mesh (zmin,zmax) "<<ZMin()<<","<<ZMax()<<endl;
  cout<<"Input::Input(): boundary (xmin) "<<Boundary(3)<<endl;
  cout<<"Input::Input(): boundary (xmax) "<<Boundary(1)<<endl;
  cout<<"Input::Input(): boundary (ymin) "<<Boundary(0)<<endl;
  cout<<"Input::Input(): boundary (ymax) "<<Boundary(2)<<endl;
  cout<<"Input::Input(): boundary (zmin) "<<Boundary(4)<<endl;
  cout<<"Input::Input(): boundary (zmax) "<<Boundary(5)<<endl;
  cout<<"Input::Input(): P (kinematics)  "<<Pk()<<endl;
  cout<<"Input::Input(): P (thermodynamics)  "<<Pt()<<endl;
  cout<<"Input::Input(): "<<NMaterials()<<" materials specified:"<<endl;
  for(int imat=0;imat<NMaterials();imat++){
  cout<<"Input::Input():      material "<<Material(imat)<<endl;
  cout<<"Input::Input():        (xmin,xmax) "<<Range(0,imat)<<","<<Range(1,imat)<<endl;
  cout<<"Input::Input():        (ymin,ymax) "<<Range(2,imat)<<","<<Range(3,imat)<<endl;
  cout<<"Input::Input():        (zmin,zmax) "<<Range(4,imat)<<","<<Range(5,imat)<<endl;
  cout<<"Input::Input():        density "<<Density(imat)<<endl;
  cout<<"Input::Input():        pressure "<<Pressure(imat)<<endl;
  cout<<"Input::Input():        velocity (x,y,z) "<<Velocity(0,imat)<<","<<Velocity(1,imat)<<","<<Velocity(2,imat)<<endl;
  }
  cout<<"Input::Input(): material "<<Ambient()<<" is set as the ambient material"<<endl;
  cout<<"Input::Input(): Done."<<endl;

  return;

}

// Accessor functions to member data

    string Input::Title(){return mTitle;} // returns problem title
    int Input::NDescriptions(){return mDescription.size();} // returns number of lines in the problem description
    string Input::Description(int i){return mDescription[i];} // returns line i of the problem description
    int Input::Nodes(int idim){return mNodes[idim];} // returns the number of nodes along direction idim
    double Input::XMin(){return mXMin;} // returns the mesh extent
    double Input::XMax(){return mXMax;} // returns the mesh extent
    double Input::YMin(){return mYMin;} // returns the mesh extent
    double Input::YMax(){return mYMax;} // returns the mesh extent
    double Input::ZMin(){return mZMin;} // returns the mesh extent
    double Input::ZMax(){return mZMax;} // returns the mesh extent
    string Input::Boundary(int iface){return mBoundaries[iface];} // returns the boundary condition of face iface
    int Input::Pk(){return mPk;} // returns the polyhedral order of the kinematic element
    int Input::Pt(){return mPt;} // returns the polyhedral order of the thermodynamic element
    int Input::NMaterials(){return mMaterial.size();} // returns the number of materials
    int Input::Material(int imat){return mMaterial[imat];} // returns material number of material imat
    int Input::Ambient(){return mAmbient;} // returns ambient material number
    double Input::Density(int imat){return mDensity[imat];} // returns density of material imat
    double Input::Pressure(int imat){return mPressure[imat];} // returns pressure in material imat
    double Input::Velocity(int idim,int imat){return mVelocity[idim][imat];} // returns velocity in direction idim of m
    double Input::Range(int idim,int imat){return mRange[idim][imat];} // return range idim of material imat
//    int Input::NDims(){return mNodes.size_of();} // return number of mesh dimensions
    int Input::NDims(){return 3;} // return number of mesh dimensions
    int Input::NCells(){return mNCells;} // return number of cells in the mesh
    bool Input::zNoh(){return mzNoh;} // returns true if mesh distortion algorithm for Noh is turned on
    bool Input::zSedov(){return mzSedov;} // returns true if mesh distortion algorithm for Sedov is turned on
    bool Input::zSaltzmann(){return mzSaltzmann;} // returns true if mesh distortion algorithm for Saltzmann is turned on

// Destructor function to release storage associated with an input class object

Input::~Input(){}

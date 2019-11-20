#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <math.h>

#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>

const double AMU = 931.4941;
const double PI = 3.1415926;
const double CUTOFF = 0.999; // cut off probability
const int ALPHA = 3; // decay mode id

// struct for one nucleus
struct Nucleus
{
  int nid;     // nucleus ID
  int spin2;   // spin*2 (doubled nuclear spin)
  int parity;  // parity 1=+, 2=-
  int isospin; // isospin 1=T<, 2=T>
  double Ex;   // Excitation energy
  double mass; // stationary mass in MeV
  double Ekin; // kinetic energy of the decaying nucleus

  // operator to compare the nuclei
  bool operator==(const Nucleus& nucl) const{
    return (nid == nucl.nid &&
	    spin2 == nucl.spin2 &&
	    parity == nucl.parity &&
	    isospin == nucl.isospin &&
	    Ex == nucl.Ex);
  }
};

// struct for one decay
struct Decay
{
  Nucleus mother;
  Nucleus daughter;
  int decaymode;
  int minj;     // minimum J transfer
  double prob;  // decay probability from the mother to the daughter
};

double calcMass( int nid );

int main(int argc, char* argv[]){
  // check the usage of the program
  if( argc !=3 ){
    std::cout << " Error for the usage of the program "
	      << std::endl;
    std::cout << "  Usage: ./evapcal <Number of events> <Excitation energy>"
	      << std::endl;
    std::cout << "Example: ./evapcal 10000 23.5" << std::endl;
    exit(-1);
  }
  // number of events to be generated
  int nevent = atoi(argv[1]);
  // 20Ne excitation energy; to open the input file
  double energy = atof(argv[2]);
  char inputfile[256], outputfile[256], rootfile[256];
  sprintf( inputfile, "./cascadeDump/DUMP_%6.3f.dat", energy);
  sprintf( outputfile, "./dat/out%6.3f.dat", energy);
  sprintf( rootfile, "./root/out%6.3f.root", energy);

  // input file from the dumped CASCADE file
  std::ifstream ifs;
  ifs.open( inputfile, std::ifstream::in);
  if( !ifs ) {
    std::cerr << " Input file opne error !" << std::endl;
    exit(-1);
  }
  // output file for the energy of alpha particle
  std::ofstream ofs;
  ofs.open( outputfile, std::ofstream::out );
  if( !ofs ) {
    std::cerr << " Output file open error !" << std::endl;
    exit(-1);
  }

  Decay oneDecay;
  std::vector<Decay> decayList;
  Nucleus mother, daughter;
  Nucleus particle;
  double decaymode;
  double prob;
  // p2pcm: momentum squared(p2) in cm system of mother
  double p2cm;
  // p2pp: momentum squared(p2) of particle in lab system
  double p2p;
  // velocity of particle in lab system
  double vp;

  // normalization factor
  double normfact;
  double dummy;

  // random number storage
  double rdecay = 0.0;
  double theta = 0.0;


  //------------------------------------------------------------------
  // for root file
  // M_ne20 - M_o16 - M_alpha
  const double qvalue = calcMass(1)-calcMass(7)-calcMass(-3);
  double mo16 = calcMass( 7 );
  double ma = calcMass( -3 );
  double ex = energy;
  double exo;

  TH1D *h1 = new TH1D("h1", "Energy spectrum of evaporating alpha",
		      181, -1.1, 35.1);
  TH1D *h2 = new TH1D("h2", "Excitation energy of {}^{16}O",
		      181, -1.1, 35.1);
  TTree *tree = new TTree("tree", "tree of evapcal");
  int totalevent = nevent;
  int evtnum=0;
  int motherid=0;
  int daughterid=0;
  int minjt=0;
  int chain=0;
  int pid = 0;
  double pEkin=0.0;
  tree->Branch("totalevent", &totalevent, "totalevent/I");
  tree->Branch("evtnum", &evtnum, "evtnum/I");
  tree->Branch("motherid", &motherid, "motherid/I");
  tree->Branch("daughterid", &daughterid, "daughterid/I");
  tree->Branch("minjt", &minjt, "minjt/I");
  tree->Branch("chain", &chain, "chain/I");
  tree->Branch("pid", &pid, "pid/I");
  tree->Branch("pEkin", &pEkin, "pEkin/D");
  tree->Branch("exo", &exo, "exo/D");
  TFile *file = new TFile( rootfile, "RECREATE");
  //------------------------------------------------------------------

  dummy = calcMass( 1 );
  // 20Ne Ex = (energy) MeV well-defined one state
  const Nucleus ne20 = {1, 0, 1, 1, energy, dummy, 0.0};

  // read the input file
  while( !ifs.eof() ){
    ifs >> mother.nid >> daughter.nid >> decaymode
	>> mother.spin2 >> daughter.spin2
	>> mother.parity >> mother.isospin
	>> daughter.parity >> daughter.isospin
	>> dummy >> dummy >> mother.Ex
	>> dummy >> dummy >> daughter.Ex
	>> prob;
    mother.mass = calcMass( mother.nid );
    daughter.mass = calcMass( daughter.nid );
    if( mother.mass<0 || daughter.mass<0 ){
      std::cout << " Unknwon nucleus in CASCADE file !"
		<< ": mother.nid = " << mother.nid << ", "
		<< ": daughter.nid = " << daughter.nid
		<< std::endl;
      exit(-1);
    }
    oneDecay.mother = mother;
    oneDecay.daughter = daughter;
    oneDecay.decaymode = decaymode;
    oneDecay.minj = (int)(fabs(daughter.spin2-mother.spin2)/2.0);
    oneDecay.prob = prob;

    // Next line is neccesary because of the specification of eof function.
    if( ifs.eof() ) break;
    decayList.push_back( oneDecay );
  }
  ifs.close();

  // set seed for the random number
  srand( (unsigned int)time(NULL) );


  // event generator loop
  for(int i=0;i<nevent;i++){
    evtnum = i+1;

    // decay chain loop
    chain = 0;
    while(1){
      chain++;
      prob = 0.0;
      normfact = 0.0;

      if( chain==1 ){
	mother = ne20;
      } else {
	mother = daughter;
      }
      
      for(std::vector<Decay>::iterator itr = decayList.begin();
	  itr != decayList.end();itr++){
	oneDecay = *itr;
	if( oneDecay.mother==mother ){
	  prob += oneDecay.prob;
	}
      }
      normfact = prob; // normalization factor
      // If the probabolity is below the cutofff,
      // the decay is not possible and skip this chain.
      if( normfact < CUTOFF ) break;
      
      // random number generator [0.0,1.0)
      rdecay = (double)rand()/RAND_MAX;
      // random decay from 20Ne
      prob = 0.0;
      for(std::vector<Decay>::iterator itr = decayList.begin();
	  itr != decayList.end();itr++){
	oneDecay = *itr;
	if( oneDecay.mother==mother ){
	  prob += oneDecay.prob/normfact;
	}
	if( rdecay < prob ){
	  break;
	}
      }

      daughter = oneDecay.daughter;
      // calculation of the kinematic of the decay
      mother.mass = calcMass( mother.nid );
      daughter.mass = calcMass( daughter.nid );
      particle.nid = - (oneDecay.decaymode);
      particle.mass = calcMass( particle.nid );
      // momentum squared in CM system of daughter nucleus
      p2cm = ((mother.Ex-daughter.Ex) + (mother.mass-daughter.mass-particle.mass))
	*(2.0*daughter.mass*particle.mass)/(daughter.mass+particle.mass);
      // random number generator [0.0,PI)
      theta = (double)rand()/RAND_MAX * PI;
      // Ekin of particle in lab system
      vp = sqrt( pow( sqrt(p2cm)/particle.mass*sin(theta), 2.0 )
		  + pow( sqrt(p2cm)/particle.mass*cos(theta)
			 - sqrt(2.0*mother.Ekin/mother.mass), 2.0 ) );
      p2p = pow(particle.mass * vp, 2.0);
      particle.Ekin = p2p / (2.0*particle.mass);
      // Ekin of daughter in lab system
      daughter.Ekin = (mother.Ex-daughter.Ex)
	+ (mother.mass-daughter.mass-particle.mass)
	- particle.Ekin;
      // kinetic  energy cutoff
      if( daughter.Ekin<0.0 ) daughter.Ekin = 0.0;
      
      // text file output
      if( oneDecay.decaymode == ALPHA ){
	ofs << std::setw(8) << i
	    << std::setw(5) << chain
	    << std::setw(5) << mother.nid
	    << std::setw(12) << std::setprecision(5)
	    << std::fixed << particle.Ekin
	    << std::endl;
	exo = ex + qvalue - particle.Ekin - ma/mo16 * particle.Ekin;
	h1->Fill( particle.Ekin );
	h2->Fill( exo );
      }
      // for root
      motherid = mother.nid;
      daughterid = daughter.nid;
      minjt = oneDecay.minj;
      pid = oneDecay.decaymode;
      pEkin = particle.Ekin;
      exo = ex + qvalue - pEkin - ma/mo16 * pEkin;
      tree->Fill();

    } // end of while loop (chain loop)

  } // end of event loop


  h1->Write();
  h2->Write();
  tree->Write();
  file->Write();
  file->Close();

  
  // close the output file
  ofs.close();
  // free the memory of vector
  std::vector<Decay>().swap(decayList);

  return 0;
}


// mass calculation. return mass in MeV from nucleus ID
double calcMass( int nid ){
  // mass excess table
  const struct{
    int nid;    // nucleus ID
    double a;   // number of nucleons
    double mex; // mass excess in MeV
  } toi[] = {
    { -1,  1.0,  8.0713 }, // neutron
    { -2,  1.0,  7.2889 }, // proton
    { -3,  4.0,  2.4249 }, // alpha
    {  1, 20.0, -7.0419 }, // 20Ne
    {  2, 19.0,  1.7520 }, // 19Ne
    {  3, 19.0, -1.4874 }, // 19F
    {  4, 18.0,  5.3176 }, // 18Ne
    {  5, 18.0,  0.8731 }, // 18F
    {  6, 18.0, -0.7828 }, // 18O
    {  7, 16.0, -4.7370 }, // 16O
    {  8, 15.0,  2.8556 }, // 15O
    {  9, 15.0,  0.1014 }, // 15N
    { 10, 12.0,  0.0    }, // 12C
    { 11, 17.0, 16.5004 }, // 17Ne
    { 12, 17.0,  1.9517 }, // 17F
    { 13, 17.0, -0.8087 }, // 17O
    { 14, 17.0, -7.8700 }, // 17N
    { 15, 14.0,  8.0077 }, // 14O
    { 16, 14.0,  2.8634 }, // 14N
    { 17, 14.0,  3.0198 }, // 14C
    { 18, 11.0, 10.6493 }, // 11C
    { 19, 11.0,  8.6677 }, // 11B
    { 20,  8.0,  4.9416 }, // 8Be
  };
  // number of the nucleus.
  // Need to edit if you add the nucleus.
  int nnucl = sizeof(toi)/sizeof(toi[0]);

  int i;
  for(i=0;i<nnucl;i++){
    if( nid == toi[i].nid ) break;
  }

  if( i >= nnucl ){
    return -1.0; // can not find the nucleus from the table.
  } else {
    return ( toi[i].a * AMU + toi[i].mex );
  }
}

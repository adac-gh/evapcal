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

const double AMU = 931.4941; // unified atomic mass unit
const double PI = 3.1415926;

// struct for one nucleus
struct Nucleus
{
  int nid;     // nucleus ID
  int spin;    // spin variable
  int parity;  // parity 1=+, 2=-
  int isospin; // isospin 1=T<, 2=T>
  double Ex;   // Excitation energy
  double mass; // stationary mass in MeV
  double Ekin; // kinetic energy for the decaying nucleus

  // operator to compare the nuclei
  bool operator==(const Nucleus& nucl) const{
    return (nid == nucl.nid &&
	    spin == nucl.spin &&
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
  double prob;  // decay probability from the mother to the daughter
};

double calcMass( int nid );

int main(int argc, char* argv[]){
  // check the usage of the program
  if( argc !=3 ){
    std::cout << " Error for the usage of the program "
	      << std::endl;
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
  Nucleus mother, daughter, gdaughter, ggdaughter;
  Nucleus particle1, particle2, particle3;
  double decaymode;
  double prob;
  double p2; // momentum squared
  // p2pcm: p2 in cm system of daughter, p2pp: p2 of particle2 in lab system
  double p2pcm, p2pp;
  // velocity of particle2 in lab system
  double vpp;
  // p2ppcm: p2 in cm system of gdaughter, p2ppp: p2 of particle3 in lab system
  double p2ppcm, p2ppp;
  // velocity of particle3 in lab system
  double vppp;

  // normalization factor
  double normfact;
  double dummy;

  // random number storage
  double rdecay1, rdecay2, rdecay3;
  double theta;


  //------------------------------------------------------------------
  // for root file
  // M_ne20 - M_o16 - M_alpha
  const double qvalue = -7.0419 - (-4.7370 + 2.4249);
  double mo16 = AMU * 16.0 - 4.7370;
  double ma = AMU * 4.0 + 2.4249;
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
  int chain=0;
  int pid = 0;
  double pEkin=0.0;
  tree->Branch("totalevent", &totalevent, "totalevent/I");
  tree->Branch("evtnum", &evtnum, "evtnum/I");
  tree->Branch("motherid", &motherid, "motherid/I");
  tree->Branch("chain", &chain, "chain/I");
  tree->Branch("pid", &pid, "pid/I");
  tree->Branch("pEkin", &pEkin, "pEkin/D");
  tree->Branch("exo", &exo, "exo/D");
  TFile *file = new TFile( rootfile, "RECREATE");
  //------------------------------------------------------------------

  dummy = calcMass( 1 );
  // 20Ne Ex = (energy) MeV well-defined one state
  const Nucleus ne20 = {1, 1, 1, 1, energy, dummy, 0.0};

  // read the input file
  while( !ifs.eof() ){
    ifs >> mother.nid >> daughter.nid >> decaymode
	>> mother.spin >> daughter.spin
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
    oneDecay.prob = prob;

    // Next line is neccesary because of the specification of eof function.
    if( ifs.eof() ) break;
    decayList.push_back( oneDecay );
  }
  ifs.close();

  // set seed for the random number
  srand( (unsigned int)time(NULL) );


  for(int i=0;i<nevent;i++){
    evtnum = i+1;
  
    // summation for the decay from 20Ne Ex = 23.5 well-defined one state
    prob = 0.0;
    normfact = 0.0;
    for(std::vector<Decay>::iterator itr = decayList.begin();
	itr != decayList.end();itr++){
      oneDecay = *itr;
      if( oneDecay.mother==ne20 ){
	prob += oneDecay.prob;
      }
    }
    normfact = prob; // normalization factor

    // random number generator [0.0,1.0)
    rdecay1 = (double)rand()/RAND_MAX;
    // random decay from 20Ne
    prob = 0.0;
    for(std::vector<Decay>::iterator itr = decayList.begin();
	itr != decayList.end();itr++){
      oneDecay = *itr;
      if( oneDecay.mother==ne20 ){
	prob += oneDecay.prob/normfact;
      }
      if( rdecay1 < prob ){
	break;
      }
    }
    mother = ne20;
    daughter = oneDecay.daughter;

    // calculation of the kinematic of the first decay
    mother.mass = calcMass( mother.nid );
    daughter.mass = calcMass( daughter.nid );
    particle1.nid = - (oneDecay.decaymode);
    particle1.mass = calcMass( particle1.nid );
    // For the first decay, 20Ne is static. So, it's simple.
    p2 = ((mother.Ex-daughter.Ex) + (mother.mass-daughter.mass-particle1.mass))
      *(2.0*daughter.mass*particle1.mass)/(daughter.mass+particle1.mass);
    daughter.Ekin = p2 / (2.0*daughter.mass);
    particle1.Ekin = p2 / (2.0*particle1.mass);
    if( oneDecay.decaymode == 3 ){
      ofs << std::setw(8) << i
	  << std::setw(5) << 1
	  << std::setw(5) << mother.nid
	  << std::setw(12) << std::setprecision(5)
	  << std::fixed << particle1.Ekin
	  << std::endl;
      exo = ex + qvalue - particle1.Ekin - ma/mo16 * particle1.Ekin;
      h1->Fill( particle1.Ekin );
      h2->Fill( exo );
    }
    // for root
    motherid = mother.nid;
    chain = 1;
    pid = oneDecay.decaymode;
    pEkin = particle1.Ekin;
    exo = ex + qvalue - pEkin - ma/mo16 * pEkin;
    tree->Fill();

    // summation for the decay from the daughter nucleus
    prob = 0.0;
    normfact = 0.0;
    for(std::vector<Decay>::iterator itr = decayList.begin();
	itr != decayList.end();itr++){
      oneDecay = *itr;
      if( oneDecay.mother==daughter ){
	prob += oneDecay.prob;
      }
    }
    normfact = prob; // normalization factor
    if( normfact < 0.999 ) {
      // Here exsits some cutoff... (due to the code CASCADE)
    } else {

      // random number generator [0.0,1.0)
      rdecay2 = (double)rand()/RAND_MAX;
      // random decay from the daughter
      prob = 0.0;
      for(std::vector<Decay>::iterator itr = decayList.begin();
	  itr != decayList.end();itr++){
	oneDecay = *itr;
	if( oneDecay.mother==daughter ){
	  prob += oneDecay.prob/normfact;
	}
	if( rdecay2 < prob ){
	  break;
	}
      }
      gdaughter = oneDecay.daughter;

      // calculation of the kinematics of the second decay
      daughter.mass = calcMass( daughter.nid );
      gdaughter.mass = calcMass( gdaughter.nid );
      particle2.nid = - (oneDecay.decaymode);
      particle2.mass = calcMass( particle2.nid );
      // momentum squared in CM system of daughter nucleus
      p2pcm = ((daughter.Ex-gdaughter.Ex) + (daughter.mass-gdaughter.mass-particle2.mass))
	*(2.0*gdaughter.mass*particle2.mass)/(gdaughter.mass+particle2.mass);
      // random number generator [0.0,PI)
      theta = (double)rand()/RAND_MAX * PI;
      // Ekin of particle2 in lab system
      vpp = sqrt( pow( sqrt(p2pcm)/particle2.mass*sin(theta), 2.0 )
		  + pow( sqrt(p2pcm)/particle2.mass*cos(theta)
			 - sqrt(2.0*daughter.Ekin/daughter.mass), 2.0 ) );
      p2pp = pow(particle2.mass * vpp, 2.0);
      particle2.Ekin = p2pp / (2.0*particle2.mass);
      // Ekin of gdaughter in lab system
      gdaughter.Ekin = (daughter.Ex-gdaughter.Ex)
	+ (daughter.mass-gdaughter.mass-particle2.mass)
	- particle2.Ekin;
      if( gdaughter.Ekin<0.0 ) gdaughter.Ekin = 0.0;
      if( oneDecay.decaymode == 3 ){
	ofs << std::setw(8) << i
	    << std::setw(5) << 2
	    << std::setw(5) << daughter.nid
	    << std::setw(12) << std::setprecision(5)
	    << std::fixed << particle2.Ekin
	    << std::endl;
	exo = ex + qvalue - particle2.Ekin - ma/mo16 * particle2.Ekin;
	h1->Fill( particle2.Ekin );
	h2->Fill( exo );
      }
      motherid = daughter.nid;
      chain = 2;
      pid = oneDecay.decaymode;
      pEkin = particle2.Ekin;
      exo = ex + qvalue - pEkin - ma/mo16 * pEkin;
      tree->Fill();

      // summation for the decay from the gdaughter nucleus
      prob = 0.0;
      normfact = 0.0;
      for(std::vector<Decay>::iterator itr = decayList.begin();
	  itr != decayList.end();itr++){
	oneDecay = *itr;
	if( oneDecay.mother==gdaughter ){
	  prob += oneDecay.prob;
	}
      }
      normfact = prob; // normalization factor
      if( normfact < 0.999 ) {
	;
      } else {
	rdecay3 = (double)rand()/RAND_MAX;
	// random decay from the gdaughter
	prob = 0.0;
	for(std::vector<Decay>::iterator itr = decayList.begin();
	    itr != decayList.end();itr++){
	  oneDecay = *itr;
	  if( oneDecay.mother==gdaughter ){
	    prob += oneDecay.prob/normfact;
	  }
	  if( rdecay3 < prob ){
	    break;
	  }
	}
	ggdaughter = oneDecay.daughter;
	gdaughter.mass = calcMass( gdaughter.nid );
	ggdaughter.mass = calcMass( ggdaughter.nid );
	particle3.nid = - (oneDecay.decaymode);
	particle3.mass = calcMass( particle3.nid );
	// momentum squared in CM system of daughter nucleus
	p2ppcm = ((gdaughter.Ex-ggdaughter.Ex) + (gdaughter.mass-ggdaughter.mass-particle3.mass))
	  *(2.0*ggdaughter.mass*particle3.mass)/(ggdaughter.mass+particle3.mass);
	// random number generator [0.0,PI)
	theta = (double)rand()/RAND_MAX * PI;
	// Ekin of particle3 in lab system
	if( gdaughter.Ekin < 0) gdaughter.Ekin = 0.0;
	vppp = sqrt( pow( sqrt(p2ppcm)/particle3.mass*sin(theta), 2.0 )
		     + pow( sqrt(p2ppcm)/particle3.mass*cos(theta)
			    - sqrt(2.0*gdaughter.Ekin/gdaughter.mass), 2.0 ) );
	p2ppp = pow(particle3.mass * vppp, 2.0);
	particle3.Ekin = p2ppp / (2.0*particle3.mass);
	// Ekin of ggdaughter in lab system
	ggdaughter.Ekin = (gdaughter.Ex-ggdaughter.Ex)
	  + (gdaughter.mass-ggdaughter.mass-particle3.mass)
	  - particle3.Ekin;
	if( ggdaughter.Ekin<0.0 ) ggdaughter.Ekin = 0.0;
	if( oneDecay.decaymode == 3 ){
	  ofs << std::setw(8) << i
	      << std::setw(5) << 3
	      << std::setw(5) << gdaughter.nid
	      << std::setw(12) << std::setprecision(5)
	      << std::fixed << particle3.Ekin
	      << std::endl;
	  exo = ex + qvalue - particle3.Ekin - ma/mo16 * particle3.Ekin;
	  h1->Fill( particle3.Ekin );
	  h2->Fill( exo );
	}
	// for root
	motherid = gdaughter.nid;
	chain = 3;
	pid = oneDecay.decaymode;
	pEkin = particle3.Ekin;
	exo = ex + qvalue - pEkin - ma/mo16 * pEkin;
	tree->Fill();
      } // end of the decision for the third decay
    } // end of the decision for the second decay
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

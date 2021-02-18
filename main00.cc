// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <fstream> // for exporting further outputs
#include <string>
using namespace Pythia8;
// Let Pythia8:: be implicit.
int main(int argc, char *argv[]) {
  int n=20; // numero de eventos
  // Begin main program.
  // Set up generation.
  Pythia pythia;
  // Declare Pythia object
  //pythia.readString("SoftQCD:inelastic = ON"); // Switch on process.
  
  std::string type = ("SoftQCD:inelastic = ON");
  
  char input=argv[1][0];
  switch (input) {
  case 's':
    type = ("SoftQCD:inelastic = ON");
    break;
  case 'h':
    type = ("HardQCD:all = ON");
    break;
  case 't':
    type = ("Top:all = ON");
    break;
  case 'w':
    type = ("WeakSingleBoson:all = ON");
  default:
    break;
  }
  

  std::cout << argv[1] << "\n";
  std::cout << type << "\n";
  
  pythia.readString(type);
  
  // you must not mix processes from the SoftQCD and HardQCD process groups, 
  // since this is likely to lead to double-counting. 
  pythia.readString("Beams:eCM = 14000."); // eV CM energy.
  std::string nse = "Next:numberShowEvent = " + std::to_string(n);
  pythia.readString(nse); // para gardar os eventos que queiramos
  
  std::string spt= "pt_.txt";
  std::string sxprod="xprod_.txt";
  std::string syprod="yprod_.txt";
  std::string szprod="zprod_.txt";
  std::string sspin="spin_.txt";
  std::string scharge="charge_.txt";
  std::string sbaryon="baryon_.txt";
  std::string set="et_.txt";
  std::ofstream outputpt(spt.insert (3,1,input)); 
  std::ofstream outputx(sxprod.insert (6,1,input)); 
  std::ofstream outputy(syprod.insert (6,1,input)); 
  std::ofstream outputz(szprod.insert (6,1,input)); 
  std::ofstream outputsp(sspin.insert (5,1,input));
  std::ofstream outputch(scharge.insert (7,1,input));
  std::ofstream outputet(set.insert (3,1,input));
  // std::ofstream outputb(sbaryon.insert (7,1,input));
  // std::ofstream outputmo("mothers.txt"); 
  // std::ofstream outputda("daughters.txt"); 
  
 
  
  
  pythia.init(); // Initialize; incoming pp beams is default.
  // Generate event(s).
  
  // int aux=0;
  
  
  for (int iEvent = 0; iEvent < n; ++iEvent) {
    pythia.next(); // Generate an(other) event. Fill event record.
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputpt << pythia.event[i].pT() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputx << pythia.event[i].xProd() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputy << pythia.event[i].yProd() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputz << pythia.event[i].zProd() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputsp << pythia.event[i].spinType() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputch << pythia.event[i].charge() << " ";
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      outputet << pythia.event[i].eT() << " ";
    }
    // for (int i = 0; i < pythia.event.size(); ++i) {
    //   outputb << pythia.event[i].isMeson() << " ";
    // }
    // for (int i = 7; i < 20; ++i) {
      
     // aux=pythia.event.daughterList(i).at(0); //haberia que facer un bucle para imprimir todas as nais
     // outputmo << aux << '\n';
     
      
    //}
   
    
    
    
  }
  (void)argc;
  return 0;
}
// End main program with error-free return.
// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <fstream> // for exporting further outputs
using namespace Pythia8;
// Let Pythia8:: be implicit.
int main() {
  int n=3; // numero de eventos
  // Begin main program.
  // Set up generation.
  Pythia pythia;
  // Declare Pythia object
  //pythia.readString("SoftQCD:inelastic = ON"); // Switch on process.
  pythia.readString("Top:all = ON");
  //pythia.readString("WeakSingleBoson:all = ON");
  
  
  // you must not mix processes from the SoftQCD and HardQCD process groups, 
  // since this is likely to lead to double-counting. 
  pythia.readString("Beams:eCM = 14000."); // eV CM energy.
  std::string nse = "Next:numberShowEvent = " + std::to_string(n);
  pythia.readString(nse); // para gardar os eventos que queiramos
  
  std::ofstream outputpt("pt.txt"); 
  std::ofstream outputx("xprod.txt"); 
  std::ofstream outputy("yprod.txt"); 
  std::ofstream outputz("zprod.txt"); 
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
    // for (int i = 7; i < 20; ++i) {
      
     // aux=pythia.event.daughterList(i).at(0); //haberia que facer un bucle para imprimir todas as nais
     // outputmo << aux << '\n';
     
      
    //}
   
    
    
    
  }
  return 0;
}
// End main program with error-free return.
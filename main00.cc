// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <fstream> // for exporting further outputs
using namespace Pythia8;
// Let Pythia8:: be implicit.
int main() {
  // Begin main program.
  // Set up generation.
  Pythia pythia;
  // Declare Pythia object
  pythia.readString("SoftQCD:inelastic = ON"); // Switch on process.
  pythia.readString("Top:all = ON");
  pythia.readString("WeakSingleBoson:all = ON");
  
  
  // you must not mix processes from the SoftQCD and HardQCD process groups, 
  // since this is likely to lead to double-counting. 
  pythia.readString("Beams:eCM = 14000."); // eV CM energy.
  pythia.readString("Next:numberShowEvent = 2"); // para gardar os eventos que queiramos
  
  std::ofstream outputpt("pt.txt"); 
  std::ofstream outputx("xprod.txt"); 
  std::ofstream outputy("yprod.txt"); 
  std::ofstream outputz("zprod.txt"); 
  
  
  pythia.init(); // Initialize; incoming pp beams is default.
  // Generate event(s).
  for (int iEvent = 0; iEvent < 2; ++iEvent) {
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
    
    
  }
  return 0;
}
// End main program with error-free return.
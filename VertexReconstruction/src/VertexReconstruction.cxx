

#include <iostream>
#include <time.h>

#include "Storage.h"
#include "Processor.h"
#include "VRAnalyzer.h"


int main(int const argc, char const * argv[]) {

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  Processor processor("FillLightEvent/pot_tree",
		      "FillLightEvent/meta_tree",
		      "FillLightEvent/event_tree",
		      {argv + 1, argv + argc});
  
  //VR
  VRAnalyzer * vrana = new VRAnalyzer("varana");
  vrana->SetVBVariables(4, 20, 50, 13, 10);
  vrana->SetProducers("pandoraNu",
		      "pandoraNu",
		      "pandoraCosmicHitRemoval",
		      "simpleFlashBeam",
		      "swtrigger");
  //vrana->RunPandora();
  vrana->RunVertexQuality();
  //vrana->RunFillTreeVariables();

  processor.AddAnalyzer(vrana);
  processor.SetOutputFileName("VertexReconstruction.root");
  processor.Run();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "\n";
  
  return 0;

}

#include <Bpp/Exceptions.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/Phylip.h>     
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Protein/JTT92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/NeighborJoining.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>
using namespace std;
using namespace bpp;

int main()
{
  try {
    const ProteicAlphabet* alphabet = new ProteicAlphabet();
    Phylip* seqReader = new Phylip(false, false);
    SiteContainer* sites = seqReader->readAlignment("Myoglobin.phy", alphabet);
    cout << "Number of sequences in alignment:      " <<  sites->getNumberOfSequences() << endl;
    cout << "Number of sites in alignment:          " <<  sites->getNumberOfSites() << endl;
    SiteContainer* completeSites = SiteContainerTools::getSitesWithoutGaps(* sites);
    cout << "Number of complete sites in alignment: " <<  completeSites->getNumberOfSites() << endl;
    delete sites;
    SubstitutionModel* model = new JTT92(alphabet);
    DiscreteDistribution* rateDist = 
        new GammaDiscreteRateDistribution(4, 0.5);
    cout << "Estimating distance matrix..." << endl; 
    DistanceEstimation distEstimation(model,
        rateDist, completeSites);
    DistanceMatrix* matrix = distEstimation.getMatrix();
    cout << "Done." << endl;
    cout << "Building tree " << endl;
    NeighborJoining nj(*matrix, false, true); // unrooted tree.
    Tree* tree = nj.getTree();
    cout << "Done." << endl;
    cout << "Estimating parameters using ML" << endl;
    DRHomogeneousTreeLikelihood* likFunction =
        new DRHomogeneousTreeLikelihood(
            *tree, *completeSites, model, rateDist, true);
    likFunction->initialize();
    // Save messages to separate files:
    ofstream* profiler = new ofstream("profile.txt", ios::out);
    ofstream* messenger = new ofstream("messages.txt", ios::out);
    OptimizationTools::optimizeNumericalParameters(likFunction,
        likFunction->getParameters(),
        0,
        5,
        0.000001, //Precision on log likelihood
        1000000,  //Max # of evaluations
        new StlOutputStreamWrapper(messenger),
        new StlOutputStreamWrapper(profiler),
        false,    //No reparametrization
        3);
    cout << "Done." << endl;
    cout << "Writing results to file 'output.dnd'." << endl;
    Newick newick;
    newick.write(likFunction->getTree(), "output.dnd");
    
    delete alphabet;
    delete seqReader;
    delete completeSites;
    delete tree;
    delete matrix;
    delete likFunction;
  }
  catch (Exception & e)
  {
    cout << "ERROR! " << e.what() << endl;
  }
}


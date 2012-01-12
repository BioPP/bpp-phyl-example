/*
 * File: ExPhylo.cpp
 * Created by: Julien Dutheil
 * Created on: Dec Wed 10 17:10 2008
 * Last modified: Jan Wed 11 09:21 2012
 * 
 * Introduction to the classes and tools for phylogenetics.
 *
 * HOW TO USE THAT FILE:
 * - General comments are written using the * * syntax.
 * - Code lines are switched off using '//'. To activate those lines, just remove the '//' characters!
 * - You're welcome to extensively modify that file!
 */

/*----------------------------------------------------------------------------------------------------*/

/*
 * We start by including what we'll need, and sort the inclusions a bit:
 */

/*
 * From the STL:
 */
#include <iostream> //to be able to output stuff in the terminal.

/*
 * We'll use the standard template library namespace:
 */
using namespace std;

/*
 * From Bpp-Seq:
 */
#include <Bpp/Seq/Alphabet.all>  /* this includes all alphabets in one shot */
#include <Bpp/Seq/Container.all> /* this includes all containers */
#include <Bpp/Seq/Io.all>        /* this includes all sequence readers and writers */

/*
 * From Bpp-Phyl:
 */
#include <Bpp/Phyl/Tree.h> /* this includes classes for tree manipulations */
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Io.all>
#include <Bpp/Phyl/Model.all> /* this includes all models */
#include <Bpp/Phyl/Distance.all>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Simulation.all>

/*
 * From Bpp-Core:
 */
#include <Bpp/Numeric/Prob.all> /* include all probability distributions */
#include <Bpp/App.all>

/*
 * All Bio++ functions are also in a namespace, so we'll use it:
 */
using namespace bpp;

/*----------------------------------------------------------------------------------------------------*/
/*
 * Now starts the real stuff...
 */


int main(int args, char ** argv)
{
  /*
   * We surround our code with a try-catch block, in case some error occurs:
   */
  try
  {
    cout << "Hello World!" << endl;

    /*
     * In this exercise we're going to work with models.
     * This mostly concerns distance-based and maximum likelihood methods.
     */

    /*
     * An important class in Bio++ Phylogenetic Library (PhylLib) is the SubstitutionModel class.
     * It describes the model of transitions between nucleotides, amino-acids or codons.
     * It is hence in tight connection with the Alphabet classes we saw before.
     *
     * Here is how it works. We will create a new instance of Kimura's 2 parameter model:
     */

    SubstitutionModel* model = new K80(&AlphabetTools::RNA_ALPHABET, 2.5);
    //SubstitutionModel* model = new JCnuc(&AlphabetTools::RNA_ALPHABET);

    /*
     * SubstitutionModel objects are quite complexe. To keep things short:
     * - They implement the Parametrizable interface, which basically means that they have Parameter
     *   objects that can be retrieved and modified.
     * - They compute the probabilities of transition between all states in the alphabet.
     * You can give a look at the documentation for classes SubstitutionModel, Parametrizable, ParameterList and Parameter.
     * Then try the following:
     */
    cout << "This model has " << model->getNumberOfStates() << " states." << endl;
    ParameterList pl = model->getParameters();
    cout << "This model has " << pl.size() << " parameters." << endl;
    pl.printParameters(cout);
    //pl.getParameter("K80.kappa").setValue(1);
    pl.printParameters(cout);
    /* Apply the parameter change: */
    //model->setParametersValues(pl);
    /* And what about that? */
    //pl.getParameter("K80.kappa").setValue(-1);

    /* INFO: The prefix "K80." is the namespace of the parameters, to avoid confusion when more complicated models are used.
     * Parameter values can also be directly accessed using:
     */
    //cout << model->getParameterValue("kappa") << endl;
    /* This syntax is then independent of the namespace used.
     */

    double t = 0.1;
    for(int i = 0; i < (int)model->getNumberOfStates(); i++)
    {
      for(int j = 0; j < (int)model->getNumberOfStates(); j++)
      {
        cout << "Probability of change from ";
        cout << model->getAlphabet()->intToChar(i) << " to ";
        cout << model->getAlphabet()->intToChar(j) << " is ";
        cout << model->Pij_t(i, j, t) << " for a branch of length ";
        cout << t << endl; 
      }
    }

    /* ----------------
     * QUESTION 1: Verify that the sum over of all Pij_t is equal to 1!
     * ----------------
     */

    /*
     * Another important set of objects are ones implementing the DiscreteDistribution interface.
     * They are also Parametrizable object, and so share a lot of methods with SubstitutionModel objects.
     * They are used here for modeling rate across site heterogeneity, but we will consider a Constant distribution for now:
     */
    DiscreteDistribution* rateDist = new ConstantDistribution(1., true);

    /*
     * We wil now use these models to build a distance matrix and a BioNJ tree:
     */
    Fasta seqReader;
    SequenceContainer* sequences = seqReader.read("SSU.fasta", &AlphabetTools::RNA_ALPHABET);
    SiteContainer* sites = new VectorSiteContainer(*sequences);
    delete sequences;
    cout << "There are " << sites->getNumberOfSites() << " positions in the alignment." << endl;
    SiteContainerTools::removeGapOnlySites(*sites);
    cout << "There are " << sites->getNumberOfSites() << " positions in the alignment that are not only made of gaps." << endl;
    SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    cout << "Computing distance matrix..." << endl;
    DistanceEstimation distanceMethod(model, rateDist, sites);
    cout << endl;
    
    /*
     * Now we retrieve the computed distances:
     */
    DistanceMatrix* distances = distanceMethod.getMatrix();

    /*
     * Now we will build a BioNJ tree:
     */
    cout << "Computing tree..." << endl;
    BioNJ bionj(*distances);
    cout << endl;
    Tree* tree = bionj.getTree();

    /*
     * And write it to a file:
     */
    Newick treeWriter;
    treeWriter.write(*tree, "SSU_BioNJ.dnd");

    /* ----------------
     * QUESTION 2: Modify the previous code in order to build a tree with a Tamura 92 model and a gamma distribution for substitution rates.
     * ----------------
     */

    /*
     * We will now use that tree to build a Maximum Likelihood tree.
     */

    NNIHomogeneousTreeLikelihood* tl = new NNIHomogeneousTreeLikelihood(*tree, *sites, model, rateDist); 
    tl->initialize();
    cout << "Log likelihood before: " << tl->getLogLikelihood() << endl;
    tl = OptimizationTools::optimizeTreeNNI(tl, tl->getParameters(), true, 100, 0.000001, 300, 1, 0, 0, false, 3);
    cout << "Log likelihood after: " << tl->getLogLikelihood() << endl;
    tl->getParameters().printParameters(cout);
    TreeTemplate<Node> mlTree(tl->getTree());
    treeWriter.write(mlTree, "SSU_ML.dnd");

    /* ----------------
     * QUESTION 3: Compare various models on this data set, for instance K80, HKY85, GTR +/- Gamma.
     * Tip: In class RandomTools (Bpp/Numeric/Random), there are tools to deal with a chi2 distribution...
     * ----------------
     */

    /*
     * Last but not least, we will now simulate data from the estimated model:
     */
    TreeTemplate<Node>* mlTreeTT = new TreeTemplate<Node>(mlTree); 
    HomogeneousSequenceSimulator* simulator = new HomogeneousSequenceSimulator(model, rateDist, mlTreeTT); 
    unsigned int numberOfSites = 500;
    SiteContainer* simSites = simulator->simulate(numberOfSites);
    Fasta seqWriter;
    seqWriter.write("Simulations.fasta", *simSites);

    /* ----------------
     * QUESTION 4: Assess some properties of the model using simulations (parametric bootstrap)
     *
     * Simulate a hundred sites using a previously fitted GTR+Gamma model,
     * then reestimate a tree and model parameters on the simulated data set.
     * Then compute the mean and variance of the estimates obtained, and compare to the real values.
     *
     * BONUS QUESTION 1: compare results between a parametric and a non-parametric boostrap.
     * BONUS QUESTION 2: also compute the consensus tree with bootstrap values (browse the documentation first!)
     * ----------------
     */

  }
  catch(Exception& e)
  {
    cout << "Bio++ exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cout << "Any other exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}


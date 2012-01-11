/*
 * File: ExTrees.cpp
 * Created by: Celine Scornavacca
 * Created on: Dec Sat 13 17:13 2008
 * Last modified: Jun Fri 05 14:21 2009 by Julien Dutheil
 *
 * Introduction to the Tree class.
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
#include <iostream> // to be able to output stuff in the terminal.

/*
 * We'll use the standard template library namespace:
 */
using namespace std;

/*
 * From the phylogenetics library:
 */
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Io.all>

using namespace bpp;

/*----------------------------------------------------------------------------------------------------*/
/*
 * Some functions we'll need:
 */
void collapseEdge(TreeTemplate<Node>& tree, Node* on)
{
  Node* temp = on->getFather();
  Node* newNode;
  unsigned i_max = (*on).getNumberOfSons();
  for (unsigned i = 0; i < i_max; i++)
  {
    newNode = on->getSon(0);   // we take always the first son...it's always a different one
    on->removeSon(on->getSon(0));
    temp->addSon(newNode);
  }
  temp->removeSon(on);
  if (temp->getNumberOfSons() == 1 && ((temp->hasFather())))
  {  // degree ==2... not so good! we have to collapse again
    collapseEdge(tree, temp);
    if ((temp->getSon(0))->hasDistanceToFather() && ( temp)->hasDistanceToFather())
      (*temp->getSon(0)).setDistanceToFather((temp->getSon(0))->getDistanceToFather() + temp->getDistanceToFather());
  }
}

void collapseOnTree(TreeTemplate<Node>& tree, float value)
{
  vector<Node*> innerNodes = tree.getInnerNodes();
  for (unsigned i = 0; i < innerNodes.size(); i++)
  {
    if (innerNodes[i]->hasFather()
        && (innerNodes[i]->hasBootstrapValue())
        && (innerNodes[i]->getBootstrapValue() < value))
    {
      collapseEdge(tree, innerNodes[i]);
    }
  }
}

/*
 * Now starts the real stuff...
 */
int main(int args, char** argv)
{
  /*
   * We surround our code with a try-catch block, in case some error occurs:
   */
  try
  {
    cout << "Hello World!" << endl;

    /*
     * In most cases, trees are represented in the newick parenthetic format.
     * Branch lengths and bootstraps are supported:
     * ex: ((A:0.1, B:0.15)90:0.2, C:0.27);
     */
    string treeDesc = "(((A:0.05,B:0.035):0.01,C:0.02):0.06,(D:0.08,E:0.08):0.11);";
    Tree* tree = TreeTemplateTools::parenthesisToTree(treeDesc);

    /*
     * Trees have branches and nodes. Terminal nodes are called leaves.
     * Trees are oriented (rooted), i.e. each node has one father node and possibly many son nodes. Leaves are
     * nodes without descendant and root is defined has the without father. Inner nodes will generally contain
     * two descendants (the tree is then called bifurcating), but multifurcating trees are also allowed with
     * this kind of description. In the rooted case, each inner node also defines a subtree. This allows to
     * work recursively on trees, which is very convenient in most cases. To deal with non-rooted trees, we
     * place an artificial root at a particular node: hence the root node appears to be trifurcated.
     */
    cout << "Number of nodes : " << tree->getNumberOfNodes() << endl;
    cout << "Number of leaves: " << tree->getNumberOfLeaves() << endl;
    vector<string> names = tree->getLeavesNames();
    for (unsigned int i = 0; i < names.size(); i++)
    {
      cout << names[i] << endl;
    }

    vector<int> ids = tree->getNodesId();
    for (unsigned int i = 0; i < ids.size(); i++)
    {
      cout << "(" << ids[i];
      if (tree->isLeaf(ids[i]))
        cout << "=" << tree->getNodeName(ids[i]);
      if (tree->isRoot(ids[i]))
        cout << "=root";
      cout << ") ";
      if (tree->hasFather(ids[i]))
        cout << "-> (" << tree->getFatherId(ids[i]) << ")";
      cout << endl;
    }

    /*
     * The ids provide a very general way to interfere with a tree.
     * Sometime, however, we need a more precise representation of the tree structure.
     * There is currently only one implementation of the Tree interface, which uses
     * special objects for representing nodes. Nodes object allow a ponter-based
     * representation of trees into memory.
     */

    /*We start creating some nodes, the fundamenta of trees*/

    /*This sets the id of node1 to 0*/
    Node* node1 = new Node(0);
    Node* node2 = new Node(1);

    /*This sets the name of node3 to species_A*/
    Node* node3 = new Node(2, "species_A");
    Node* node4 = new Node(3, "species_B");
    Node* node5 = new Node(4, "species_C");

    /*Let's check!*/

    cout << "The id of node1 is " << node1->getId() << endl;
    cout << "The name of node3 is " << node3->getName() << endl;

    /* ----------------
     * QUESTION 1: - Try to print the name of node2. What's the problem?
     *             - Change the names of node3, node4, node5 to true species
     * ----------------
     */

    /*We create parental relationships between these nodes*/

    node1->addSon(node2);
    node1->addSon(node3);
    node2->addSon(node4);
    node2->addSon(node5);

    /*
     * Then we create a tree from these nodes.
     * In order to do so, we use the TreeTemplate class,
     * which implements the Tree interface.
     */

    TreeTemplate<Node>* tree1 = new TreeTemplate<Node>(node1);
    tree1->setName("Tree1");
    cout << "The tree " << tree1->getName() << " has " << tree1->getNumberOfNodes() << " nodes" << endl;
    cout << "The tree " << tree1->getName() << " has " << tree1->getNumberOfLeaves() << " leaves" << endl;

    /*
     * NB: the parenthesisToTree method that we used before actually creates a TreeTemplate<Node> object !
     */


    /*
     * In most cases, trees are created from a file and not from scratch.
     * We hence need to introduce a new kind of objects, called newick.
     * These objects have a method 'read' that creates a tree from a file and
     * a method 'write' that prints a tree on a specified output.
     */

    Newick treeIO(true);

    cout << "The tree " << tree1->getName() << " is: ";
    treeIO.write(*tree1, cout);

    /* ----------------
     * QUESTION 2: - Uncomment the following lines and print the tree again. Has the tree a new leaf? Why?
     * ----------------
     */

    Node* node6 = new Node(5, "species_D");
    node6->setFather(node2);
    treeIO.write(*tree1, cout);
    delete tree1;

    /* ----------------
     * QUESTION 3: - Try to create a new tree tree2 from the file tree.txt. Hint: use the newick object
     *               You can visualise this tree using NJplot, TreeView, dendroscope....
     * ----------------
     */
    Tree* tmp = treeIO.read("tree.txt");
    TreeTemplate<Node>* tree2 = new TreeTemplate<Node>(*tmp);
    delete tmp;

    /* We can find a node in a tree looking for its name or id: they need to be unique!*/
    cout << "The name of the taxon whose name is Zea is: " << (tree2->getNode("Zea"))->getName() << "... bravo!\n";

    /* ----------------
     * QUESTION 4: - Reroot tree2 on the taxon whose name is Brachyp. Hint: you need to define a new outgroup
     *               Print tree2 to check
     * ----------------
     */

    /* ----------------
     * QUESTION 5: - Create a new tree tree3 from the file tree2.txt and use the Grafen's method to initialize branch lengths.
     * ----------------
     */

    /*
     * We can combine trees into a larger one. We add tree3 as subtree son of the father of the taxon whose name is Zea.
     * We need to reset the id of all nodes because, merging the two trees, we can have an id more than once!
     */

    TreeTemplate<Node>* tree3 = treeIO.read("tree2.txt");
    Node* temp = tree2->getNode("Zea")->getFather();
    temp->addSon(tree3->getRootNode());
    tree2->resetNodesId();
    treeIO.write(*tree2, cout);

    /*We can also make a copy of a tree*/

    TreeTemplate<Node>* tree4 = tree2->clone();

    /* ----------------
     * QUESTION 6: - Try to  modify tree4, change taxon names or branch lenghts.... Use the documentation of the class
     *               TreeTemplate and TreeTemplateTools to find some ideas. Have fun!
     *               Then print tree2 and tree4. What do you notice?
     * ----------------
     */

    /*
     * The function collapseOnTree is an example of how to manipulate trees in Bio++.
     * It collapses all branches with a boostrap value < to a given threshold
     */

    float threshold = 10;
    collapseOnTree(*tree2, threshold);
    treeIO.write(*tree2, cout);
    delete tree2;

    /* ----------------
     * QUESTION 7: - Try to understand how collapseOnTree works and define a new function that collapses
     *               all branches that are shorter than a threshold.
     * ----------------
     */

    /*Let's read more than a tree in only one line!*/

    vector<Tree*> vecTr;
    treeIO.read("Orthomam_v1_787trees.trees", vecTr);

    /*
     * With bio++ we can easly compute (strict and majority) consensus trees and MRP supertrees
     */

    TreeTemplate<Node>* strictConsensusTree = TreeTools::strictConsensus(vecTr, false);
    TreeTemplate<Node>* majorityConsensusTree = TreeTools::majorityConsensus(vecTr, false);

    /*
     * Have a look to the documentation: why I defined the result of TreeTools::MRP as an object Tree and
     * not as an object TreeTemplate<Node>?
     */

    Tree* MRPTree = TreeTools::MRP(vecTr);

    /* ----------------
     * QUESTION 8: Compute the Robinson and Foulds Distance of majorityConsensusTree and MRPTree
     * ----------------
     */

    delete strictConsensusTree;
    delete majorityConsensusTree;
    delete MRPTree;

    for (unsigned int i = 0; i < vecTr.size(); i++)
    {
      delete vecTr[i];
    }

    /*
     * This is only a tiny tiny introduction to the data structures for phylogenetics that Bio++ offers. Among other things,
     * an interesting feature is the possibility to define your personal data structure, adding new
     * attributes to the standard class Node using the NodeTemplate<NodeInfos> Class Template
     */
  }
  catch (Exception& e)
  {
    cout << "Bio++ exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }
  catch (exception& e)
  {
    cout << "Any other exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}


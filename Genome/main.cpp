#include <iostream>
#include <fstream>
#include "provided.h"
#include "Trie.h"
using namespace std;

int main()
{
    Trie<int>* t = new Trie<int>();
    assert(t->find("Tree", true).size() == 0);
    assert(t->find("Tree", false).size() == 0);
    assert(t->find("Trie", true).size() == 0);
    assert(t->find("Trie", false).size() == 0);
    t->insert("Tree", 5);
    assert(t->find("Tree", true).size() == 1);
    assert(t->find("Tree", false).size() == 1);
    assert(t->find("Trie", true).size() == 0);
    assert(t->find("Trie", false).size() == 1);
    t->insert("Tree", 5);
    assert(t->find("Tree", true).size() == 2);
    assert(t->find("Tree", false).size() == 2);
    assert(t->find("Trie", true).size() == 0);
    assert(t->find("Trie", false).size() == 2);
    t->insert("Trie", 6);
    assert(t->find("Tree", true).size() == 2);
    assert(t->find("Tree", false).size() == 3);
    assert(t->find("Trie", true).size() == 1);
    assert(t->find("Trie", false).size() == 3);
    t->insert("Trees", 7);
    assert(t->find("Tree", true).size() == 2);
    assert(t->find("Tree", false).size() == 3);
    assert(t->find("Trie", true).size() == 1);
    assert(t->find("Trie", false).size() == 3);
    cout << "All Trie tests passed!" << endl;
    
    string filename = "/Users/Richard/Downloads/Genomics\ Data/Halorubrum_chaoviator.txt";
    ifstream strm(filename);
    if (!strm)
    {
        cout << "Cannot open " << filename << endl;
        return 1;
    }
    vector<Genome> vg;
    bool success = Genome::load(strm, vg);
    
    if (success)
    {
        cout << "Loaded " << vg.size() << " genomes successfully:" << endl;
        for (int k = 0; k != vg.size(); k++)
            cout << vg[k].name() << ": " << vg[k].length() << " bases" << endl;
    }
    else
        cout << "Error loading genome data" << endl;
    
}

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
    t->reset();
    assert(t->find("Tree", true).size() == 0);
    assert(t->find("Tree", false).size() == 0);
    assert(t->find("Trie", true).size() == 0);
    assert(t->find("Trie", false).size() == 0);
    
    //
    // Test Trie copying, and assignment
    //
    
    cout << "All Trie tests passed!" << endl;
    
    //string filename = "/Users/Richard/Downloads/Genomics\ Data/Halorubrum_chaoviator.txt";
    string filename = "/Users/Richard/Documents/CS32/Genome/Genome/Text.txt";
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
    
    GenomeMatcher matcher(4);
    assert(matcher.minimumSearchLength() == 4);
    matcher.addGenome(vg[0]);
    matcher.addGenome(vg[1]);
    matcher.addGenome(vg[2]);
    vector<DNAMatch> matches;
    bool result;
    
    result = matcher.findGenomesWithThisDNA("GAAG", 4, true, matches);
    assert(result && matches.size() == 3);
    assert(matches[0].length == 4 && matches[0].position == 60);
    assert(matches[1].length == 4 && matches[1].position == 54);
    assert(matches[2].length == 4 && matches[2].position == 29);
    
    result = matcher.findGenomesWithThisDNA("GAATAC", 4, true, matches);
    assert(result && matches.size() == 2);
    assert(matches[0].length == 5 && matches[0].position == 22);
    assert(matches[1].length == 5 && matches[1].position == 48);
    
    result = matcher.findGenomesWithThisDNA("GAATAC", 6, true, matches);
    assert(!result && matches.size() == 0);
    
    result = matcher.findGenomesWithThisDNA("GAATAC", 6, false, matches);
    assert(result && matches.size() == 2);
    assert(matches[0].length == 6 && matches[0].position == 22);
    assert(matches[1].length == 6 && matches[1].position == 48);
    
    result = matcher.findGenomesWithThisDNA("GTATAT", 6, false, matches);
    assert(result && matches.size() == 2);
    assert(matches[0].length == 6 && matches[0].position == 22);
    assert(matches[1].length == 6 && matches[1].position == 48);
    
    result = matcher.findGenomesWithThisDNA("GAATACG", 6, false, matches);
    assert(result && matches.size() == 2);
    assert(matches[0].length == 6 && matches[0].position == 22);
    assert(matches[1].length == 7 && matches[1].position == 48);
    
    result = matcher.findGenomesWithThisDNA("GAAGGGTT", 5, false, matches);
    assert(result && matches.size() == 3);
    assert(matches[0].length == 8 && matches[0].position == 60);
    assert(matches[1].length == 5 && matches[1].position == 54);
    assert(matches[2].length == 7 && matches[2].position == 35);
    
    result = matcher.findGenomesWithThisDNA("GAAGGGTT", 6, false, matches);
    assert(result && matches.size() == 2);
    assert(matches[0].length == 8 && matches[0].position == 60);
    assert(matches[1].length == 7 && matches[1].position == 35);
    
    result = matcher.findGenomesWithThisDNA("ACGTGCGAGACTTAGAGCC", 12, false, matches);
    assert(result && matches.size() == 1);
    assert(matches[0].length == 19 && matches[0].position == 28);
    
    result = matcher.findGenomesWithThisDNA("ACGTGCGAGACTTAGAGCG", 12, false, matches);
    assert(result && matches.size() == 1);
    assert(matches[0].length == 19 && matches[0].position == 28);
    
    result = matcher.findGenomesWithThisDNA("GAAG", 3, true, matches);
    assert(!result && matches.size() == 0);
    
    result = matcher.findGenomesWithThisDNA("GAAG", 5, true, matches);
    assert(!result && matches.size() == 0);
    
    //cout << (result ? "True" : "False") << endl;
    //for (int i = 0; i < matches.size(); i++)
    //    cout << matches[i].genomeName << " of length " << matches[i].length << " at position " << matches[i].position << endl;
    
    Genome example("Example", "ACGTGCGAGACTTAGAGCC");
    vector<GenomeMatch> results;
    bool s;
    s = matcher.findRelatedGenomes(example, 12, true, 0.99999, results);
    assert(s && results.size() == 1 && results[0].genomeName == "Genome 2");
    s = matcher.findRelatedGenomes(example, 1, true, 0.99999, results);
    assert(!s);
    s = matcher.findRelatedGenomes(example, 5, false, 0, results);
    assert(s && results.size() == 3);
    //for (int i = 0; i < results.size(); i++)
    //    cout << results[i].percentMatch << endl;
    
    cout << "All Genome Matcher tests passed!" << endl;
}

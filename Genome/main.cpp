#include <iostream>
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
    cout << "All tests passed!" << endl;
}

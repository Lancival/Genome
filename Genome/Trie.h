#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    class TrieNode
    {
    public:
        TrieNode();
        ~TrieNode();
        void insert(const std::string& key, const ValueType& value, const int& i);
        void find(const std::string& key, bool exactMatchOnly, const int& i, std::vector<ValueType>& answer) const;
    private:
        std::vector<ValueType> values;
        std::vector<char> labels;
        std::vector<TrieNode*> children;
    };
    TrieNode* root;
};

template<typename ValueType>
Trie<ValueType>::TrieNode::TrieNode() : values(), labels(), children() {}

template<typename ValueType>
Trie<ValueType>::TrieNode::~TrieNode()
{
    for (int i = 0; i < children.size(); i++)
        delete children[i];
}

template<typename ValueType>
void Trie<ValueType>::TrieNode::insert(const std::string& key, const ValueType& value, const int& i)
{
    if (i == key.size())
        return values.push_back(value);
    for (int l = 0; l < labels.size(); l++)
        if (labels[l] == key[i])
            return children[l]->insert(key, value, i+1);
    labels.push_back(key[i]);
    children.push_back(new TrieNode());
    children[children.size() - 1]->insert(key, value, i+1);
}

template<typename ValueType>
void Trie<ValueType>::TrieNode::find(const std::string& key, bool exactMatchOnly, const int &i, std::vector<ValueType>& answer) const
{
    if (i == key.size())
        for (int i = 0; i < values.size(); i++)
            answer.push_back(values[i]);
    else
    {
        for (int l = 0; l < labels.size(); l++)
        {
            if (labels[l] == key[i])
                children[l]->find(key, exactMatchOnly, i+1, answer);
            else if (!exactMatchOnly && i != 0)
                children[l]->find(key, true, i+1, answer);
        }
    }
}

template<typename ValueType>
Trie<ValueType>::Trie()
{
    root = new TrieNode();
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
    delete root;
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    delete root;
    root = new TrieNode();
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    root->insert(key, value, 0);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    std::vector<ValueType> answer;
    root->find(key, exactMatchOnly, 0, answer);
    return answer;
}

#endif // TRIE_INCLUDED

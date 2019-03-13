#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength;
    vector<Genome> genomes;
    Trie<pair<int,int>> sequences;
    unordered_map<string,int> names;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength) : genomes(), sequences(), names()
{
    m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    genomes.push_back(genome);
    names.insert({genome.name(), genomes.size() - 1});
    for (int i = 0; i <= genome.length() - m_minSearchLength; i++)
    {
        string fragment;
        genome.extract(i, m_minSearchLength, fragment);
        sequences.insert(fragment, {genomes.size() - 1, i});
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    matches.clear();
    if (fragment.length() < minimumLength || minimumLength < m_minSearchLength)
        return false;
    vector<pair<int, int>> potentialMatches = sequences.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);
    unordered_map<int, int> length;
    unordered_map<int, int> position;
    for (int i = 0; i < potentialMatches.size(); i++)
    {
        bool exactMatch = true;
        int g = potentialMatches[i].first;
        int p = potentialMatches[i].second;
        for (int j = 0; j < fragment.length(); j++)
        {
            string base;
            genomes[g].extract(p + j, 1, base);
            if (base[0] != fragment[j] && !exactMatchOnly && exactMatch)
                exactMatch = false;
            else if (base[0] != fragment[j])
                break;
            if (j+1 >= minimumLength && (length[g] < j+1 || (length[g] == j+1 && position[g] > p)))
            {
                length[g] = j+1;
                position[g] = p;
            }
        }
    }
    for (int g = 0; g < genomes.size(); g++)
        if (length.find(g) != length.end())
            matches.push_back({genomes[g].name(), length[g], position[g]});
    return (matches.size() == 0) ? false : true;
}

bool compareGenomeMatch(GenomeMatch a, GenomeMatch b)
{
    return (a.percentMatch == b.percentMatch) ? (a.genomeName < b.genomeName) : (a.percentMatch > b.percentMatch);
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength < m_minSearchLength)
        return false;
    results.clear();
    unordered_map<int, int> matches;
    for (int i = 0; i < query.length()/fragmentMatchLength; i++)
    {
        string fragment;
        query.extract(i*fragmentMatchLength, fragmentMatchLength, fragment);
        vector<DNAMatch> m;
        findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, m);
        for (int j = 0; j < m.size(); j++)
            matches[names.at(m[j].genomeName)]++;
    }
    for (unordered_map<int, int>::iterator it = matches.begin() ; it != matches.end(); it++)
    {
        double percentMatch = ((double)(it->second))/(query.length()/fragmentMatchLength)*100;
        if (percentMatch > matchPercentThreshold)
            results.push_back({genomes[it->first].name(), percentMatch});
    }
    sort(results.begin(), results.end(), compareGenomeMatch);
    return (results.size() == 0) ? false : true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}

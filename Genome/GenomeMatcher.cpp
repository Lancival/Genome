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
    Trie<string> sequences;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength) : genomes(), sequences()
{
    m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    genomes.push_back(genome);
    for (int i = 0; i <= genome.length() - m_minSearchLength; i++)
    {
        string fragment;
        genome.extract(i, m_minSearchLength, fragment);
        sequences.insert(fragment, to_string(genomes.size()-1) + "," + to_string(i));
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
    vector<string> potentialMatches = sequences.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);
    unordered_map<int, int> length;
    unordered_map<int, int> position;
    for (int i = 0; i < potentialMatches.size(); i++)
    {
        bool exactMatch = true;
        size_t delim = potentialMatches[i].find(',');
        int g = stoi(potentialMatches[i].substr(0, delim));
        int p = stoi(potentialMatches[i].substr(delim+1, potentialMatches[i].size()));
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
    vector<int> matches;
    for (int i = 0; i < genomes.size(); i++)
        matches.push_back(0);
    for (int i = 0; i < query.length()/fragmentMatchLength; i++)
    {
        string fragment;
        query.extract(i*fragmentMatchLength, fragmentMatchLength, fragment);
        vector<DNAMatch> m;
        findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, m);
        for (int j = 0, g = 0; j < m.size(); j++)
        {
            while (m[j].genomeName != genomes[g].name())
                g++;
            matches[g]++;
        }
    }
    for (int g = 0; g < genomes.size(); g++)
    {
        double percentMatch = ((double)(matches[g]))/(query.length()/fragmentMatchLength)*100;
        if (percentMatch > matchPercentThreshold)
        {
            GenomeMatch gm;
            gm.genomeName = genomes[g].name();
            gm.percentMatch = percentMatch;
            results.push_back(gm);
        }
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

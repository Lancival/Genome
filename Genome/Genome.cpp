#include "provided.h"
#include <string>
#include <vector>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
    size_t m_size;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
    m_size = m_sequence.size();
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    genomes.clear();
    char c;
    genomeSource.get(c);
    if (c != '>')
        return false;
    while (!genomeSource.eof())
    {
        string name;
        string sequence;
        getline(genomeSource, name);
        if (name.size() == 0)
            return false;
        while (!genomeSource.eof())
        {
            genomeSource.get(c);
            c = toupper(c);
            if (c == '\n' && genomeSource.peek() != '\n')
                continue;
            else if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
                sequence += c;
            else if (c == '>' || c == EOF)
                break;
            else
                return false;
        }
        if (sequence.size() == 0)
            return false;
        genomes.push_back(Genome(name, sequence));
    }
    return true;
}

int GenomeImpl::length() const
{
    return (int)(m_size);
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if (position < 0 || position > m_size - length)
        return false;
    fragment = m_sequence.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}

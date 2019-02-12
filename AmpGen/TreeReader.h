#ifndef AMPGEN_TREEREADER_H
#define AMPGEN_TREEREADER_H 1

#include "AmpGen/MsgService.h"
#include "TLeaf.h"
#include "TTree.h"
#include <vector>

namespace AmpGen
{
  template <class OutputType>
    class TreeReader 
    {
      private:
        struct IReadBranch {
          std::string name;
          IReadBranch( const std::string& name = "" ) : name( name ) {}
          virtual void* address() const        = 0;
          virtual void transfer()              = 0;
          virtual ~IReadBranch()               = default;
        };

        template <class InputType> struct ReadBranch : public IReadBranch {
          InputType thing;
          OutputType* output;
          void* address() const override { return (void*)&thing; }
          ReadBranch( const std::string& name, OutputType* outputBranch ) : IReadBranch( name ), output( outputBranch ) {}
          void transfer() override { *output = thing; }
        };
  
       template <class InputType> struct ReinterpretBranch : public IReadBranch {
          InputType thing;
          OutputType* output;
          void* address() const override { return (void*)&thing; }
          ReinterpretBranch( const std::string& name, OutputType* outputBranch ) : IReadBranch( name ), output( outputBranch ) {}
          void transfer() override { *output = reinterpret_cast<OutputType>(thing); }
        };

        struct Branch {
          OutputType* value;
          Branch() : value( new OutputType() ) {}
          ~Branch() { delete value; }
          operator OutputType() const { return *value; }
          operator OutputType&() { return *value; }
          OutputType* operator&() { return value; }
        };

        struct Iterator {
          size_t m_position;
          TreeReader* m_parent;
          Iterator( const size_t& pos, TreeReader* parent ) : m_position( pos ), m_parent( parent ) {}
          Iterator& operator++()
          {
            m_position++;
            m_parent->getEntry( m_position );
            return *this;
          }
          bool operator==( const Iterator& rhs ) const { return m_position == rhs.m_position; }
          bool operator!=( const Iterator& rhs ) const { return m_position != rhs.m_position; }
          size_t operator*() const { return m_position; }
        };
        TTree* tree;
        bool ready;
        std::vector<IReadBranch*> branches;

      public: 
        TreeReader( TTree* tree ) : tree( tree ) {}
        void setBranch( const std::string& name, OutputType* ptr )
        {
          IReadBranch* new_branch                    = nullptr;
          TLeaf* leaf = tree->GetLeaf( name.c_str() );
          if( leaf == nullptr ){
            ERROR( "Leaf: " << name << " not found");
            return; 
          }
          std::string branchType                     = leaf->GetTypeName();
          if( branchType == "Double_t" ) new_branch = new ReadBranch<Double_t>( name, ptr );
          if( branchType == "Float_t" )  new_branch = new ReadBranch<Float_t>( name, ptr );
          if( branchType == "Bool_t" )   new_branch = new ReadBranch<Bool_t>( name, ptr );
          if( branchType == "Int_t" )    new_branch = new ReadBranch<Int_t>( name, ptr );
          if( branchType == "UInt_t" )   new_branch = new ReadBranch<UInt_t>( name, ptr );
          if( branchType == "ULong64_t") new_branch = new ReadBranch<ULong64_t>( name, ptr );
          if( new_branch == nullptr ){
            ERROR( "Branch type:" << branchType <<  " not recognised" );
            return;   
          }
          ready = false;
          branches.push_back( new_branch );
        }
        Branch bind( const std::string& name )
        {
          Branch rt;
          setBranch( name, &rt );
          return rt;
        }
        void getEntry( const unsigned int& entry )
        {
          if ( !ready ) prepare();
          tree->GetEntry( entry );
          for ( auto& branch : branches ) branch->transfer();
        }
        void prepare()
        {
          for ( auto& branch : branches ) {
            tree->SetBranchStatus( branch->name.c_str(), "1" );
            tree->SetBranchAddress( branch->name.c_str(), branch->address() );
          }
          ready = true;
        }
        ~TreeReader()
        {
          for ( auto& branch : branches ) delete branch;
        }
        Iterator begin() { return Iterator( 0, this ); }
        Iterator end() { return Iterator( tree->GetEntries(), this ); }
    };
} // namespace AmpGen

#endif

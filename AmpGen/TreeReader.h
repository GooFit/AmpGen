#ifndef AMPGEN_TREEREADER_H
#define AMPGEN_TREEREADER_H 1

#include "AmpGen/MsgService.h"
#include "TLeaf.h"
#include "TTree.h"
#include <vector>
#include "AmpGen/MetaUtils.h"
namespace AmpGen
{
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

      template <class InputType, class OutputType> struct ReadBranch : public IReadBranch 
      {
        InputType thing;
        OutputType* output;
        void* address() const override { return (void*)&thing; }
        ReadBranch( const std::string& name, OutputType* outputBranch ) : IReadBranch( name ), output( outputBranch ) {}
        void transfer() override { *output = thing; }
      };

      template <class InputType, class OutputType> struct ReinterpretBranch : public IReadBranch 
      {
        InputType thing;
        OutputType* output;
        void* address() const override { return (void*)&thing; }
        ReinterpretBranch( const std::string& name, OutputType* outputBranch ) : IReadBranch( name ), output( outputBranch ) {}
        void transfer() override { *output = reinterpret_cast<OutputType>(thing); }
      };
      struct Iterator {
        size_t      m_position;
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
        unsigned pos() const { return m_position; }
      };
      TTree* m_tree                        = {nullptr};
      bool   m_ready                       = {false};
      std::vector<IReadBranch*> m_branches = {};
      std::vector<size_t> m_entryList    = {};
    public: 
      template <class OutputType> struct Proxy 
      {
        OutputType* data;
        Proxy() : data( new OutputType() ) {} 
        operator OutputType() const { return *data; } 
        ~Proxy(){ delete data ; }
      };
      explicit TreeReader( TTree* tree );
      template <typename OutputType> Proxy<OutputType> make_proxy( const std::string& name )
      {
        Proxy<OutputType> rt;
        setBranch( name, rt.data );
        return rt;   
      }
      template <typename OutputType> void setBranch( const std::string& name, OutputType& thing ) {
       setBranch(name,&thing) ; } 
      template <typename OutputType> void setBranch( const std::string& name, OutputType* ptr )
      {
        IReadBranch* new_branch = nullptr;
        TLeaf* leaf = m_tree->GetLeaf( name.c_str() );
        if( leaf == nullptr ){
          ERROR( "Leaf: " << name << " not found");
          return; 
        }
        std::string branchType                     = leaf->GetTypeName();
        if( branchType == "Double_t" ) new_branch = new ReadBranch<Double_t , OutputType>( name, ptr );
        if( branchType == "Float_t" )  new_branch = new ReadBranch<Float_t  , OutputType>( name, ptr );
        if( branchType == "Bool_t" )   new_branch = new ReadBranch<Bool_t   , OutputType>( name, ptr );
        if( branchType == "Int_t" )    new_branch = new ReadBranch<Int_t    , OutputType>( name, ptr );
        if( branchType == "UInt_t" )   new_branch = new ReadBranch<UInt_t   , OutputType>( name, ptr );
        if( branchType == "ULong64_t") new_branch = new ReadBranch<ULong64_t, OutputType>( name, ptr );
        if( new_branch == nullptr ){
          ERROR( "Branch type:" << branchType <<  " not recognised" );
          return;   
        }
        DEBUG("Making branch with properties: [name = " << name  << ", input type = "  << branchType << " output type = " << type_string<OutputType>() << "]" );
        m_ready = false;
        m_branches.push_back( new_branch );
      }
      void setEntryList( const std::vector<size_t>& entryList );
      void unsetEntryList(); 
      void getEntry( const unsigned int& entry );
      void prepare();
      size_t nEntries() const;
      ~TreeReader();
      Iterator begin();
      Iterator end();
  };
//  ENABLE_DEBUG(TreeReader);
} // namespace AmpGen

#endif

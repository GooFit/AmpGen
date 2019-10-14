#ifndef QCGenerator
#define QCGenerator

namespace AmpGen{

class QCGenerator : public Generator{
    private : 
        EventType m_tagType;
    public : 
        void FillEventList(int nEvents, EventList list1, EventList list2);

};

}
#endif

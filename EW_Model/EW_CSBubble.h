#ifndef EW_CSBUBBLE_H
#define EW_CSBUBBLE_H

#include <cmath>
#include "EW_BubbleNucl.h"
#include "EW_Winding.h"

namespace EW_BubbleNucleation {
    using namespace Electroweak;

    template<int DIM>
    class CSBubble : public BubbleNucleation<DIM> {
    public:
        CSBubble(const Lattice<DIM>& lat, 
                const Parallel2D& parallel, 
                const std::string& id);
        ~CSBubble() = default;

        void OneBubble_Transformed(const int winding) const;
    protected:

    };

    template<int DIM>
    CSBubble<DIM>::CSBubble(
        const Lattice<DIM>& lat, 
        const Parallel2D& parallel, 
        const std::string& id)
        :
        BubbleNucleation<DIM>(lat, parallel, id)
    {}

    template<int DIM>
    void CSBubble<DIM>::OneBubble_Transformed(const int winding){
        
    }
}


#endif
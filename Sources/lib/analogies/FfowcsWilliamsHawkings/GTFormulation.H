#ifndef GTFormulation_H
#define GTFormulation_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class GTFormulation : public fwhFormulation
{

private:

    //-
    List<List<tdPrimitive<vector> > > Q_;

    //-
    List<List<tdPrimitive<vector> > > L_;
    
private:
    
    //- Disable default constructor
    GTFormulation();
    
    //- Disable copy constructor
    GTFormulation(const GTFormulation&);

public:
    
    //- Constructs a new GTFormulation object using
    //the reference to FfowcsWilliamsHawkings object
    GTFormulation(const FfowcsWilliamsHawkings& fwh);
    
    //- Destructor for an GTFormulation object
    ~GTFormulation();

    //- Calculates acoustic data used to compute acoustic pressure
    // for the specified observer by employing GT formulation of FWH Analogy
    void calculateAcousticData
    (
        const vectorField& Sf,
        const vectorField& uS,
        const scalarField& rhoS,
        const scalarField& pS,
        label iObs,
        label iSurf,
        scalar ct
    ) override;
    
    //-
    void update();
    
    //-
    void initialize();
};

}

} // namespace Foam

#endif

// ************************************************************************* //

//
//END OF FILE
//
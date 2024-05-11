#ifndef ReadFromFile_H
#define ReadFromFile_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "AcousticAnalogy.H"
#include "IFstream.H"

namespace Foam
{

namespace functionObjects
{

class ReadFromFile : public fwhFormulation
{

private:

    //-
    autoPtr<IFstream> qdsInput_;

    //-
    autoPtr<IFstream> fdsInput_;

    //-
    const label n_blocks_to_read_;
    
private:
    
    //- Disable default constructor
    ReadFromFile();
    
    //- Disable copy constructor
    ReadFromFile(const ReadFromFile&);

public:
    
    //-
    ReadFromFile(const FfowcsWilliamsHawkings& fwh);
    
    //-
    ~ReadFromFile();
    
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
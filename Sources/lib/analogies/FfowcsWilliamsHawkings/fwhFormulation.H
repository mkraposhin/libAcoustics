#ifndef fwhFormulation_H
#define fwhFormulation_H

#include "label.H"
#include "scalar.H"
#include "DynamicList.H"
#include "runTimeSelectionTables.H"
#include "AcousticAnalogy.H"
#include "OFstream.H"

namespace Foam
{

namespace functionObjects
{

class FfowcsWilliamsHawkings;

#define forAllSurfaces(arr, iObs, iSurf) \
forAll(arr, iObs)                        \
forAll(arr[iObs], iSurf)

#define forAllFaces(arr, iObs, iSurf, iFace) \
forAllSurfaces(arr, iObs, iSurf)             \
forAll(arr[iObs][iSurf], iFace)

class fwhFormulation
{

protected:

    //types for internal usage
    typedef Pair<DynamicList<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    
    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

protected:
    
    //-
    const FfowcsWilliamsHawkings& fwh_;
    
    //-
    label fwhProbeI_;

    //-
    List<List<label>> dumpProbeI_;

    //-
    surfaceTimeData qds_;
    
    //-
    surfaceTimeData fds_;

    //-
    tdPrimitive<scalar> intDotQdS_;
    
    //-
    tdPrimitive<scalar> intFdS_;
    
    //- Keys for identifying observer time
    List<List<List<scalar> > > tobs_;
    
    //- Distances from observers to surfaces elements
    List<List<List<vector> > > robs_;
    
    //- Mag of distances from observers to surfaces elements
    List<List<List<scalar> > > magrobs_;
    
    //-
    List<tdPrimitive<vector> > ni_;
    
    //-
    List<List<scalar> > nl_;

    //-
    List<scalar> rMax_;
    
    //-
    List<scalar> tauMax_;
    
    //-
    List<scalar> rMin_;
    
    //-
    List<scalar> tauMin_;

    //-
    word qdsName_;

    //-
    word fdsName_;

    //-
    autoPtr<OFstream> qdsDat_;

    //-
    autoPtr<OFstream> fdsDat_;
    

protected:
    
    //- Initializes the object
    void initialize();

    //- Returns the reference to a qds file
    OFstream& qds_out(label iSurf, label iObs);

    //- Returns the reference to an fds file
    OFstream& fds_out(label iSurf, label iObs);

    //- Returns the name of qds data file
    const word& qdsName() const;

    //- Returns the name of fds data file
    const word& fdsName() const;

private:
    //- forbid default constructor
    fwhFormulation();

    //- forbid copy constructor
    fwhFormulation(const fwhFormulation&);

public:
    
    //- constructs using given analogy object
    fwhFormulation(const FfowcsWilliamsHawkings& analogy);
    
    //-
    virtual ~fwhFormulation();
    
    //-
    scalar valueAt (const surfaceTimeData& data, label iObs, label iSurf, label iFace, scalar tau);
    
    //-
    label findExpiredIndex (const pointTimeData& timeData, scalar expiredTime);
    
    //-
    pointTimeData getNewPointData (const pointTimeData& timeData, label expiredIndex);
    
    //-
    surfaceTimeData getQdsData();
    
    //-
    surfaceTimeData getFdsData();
    
    //- Calculates acoustic pressure at the specified time
    scalar observerAcousticPressure
    (
        const vectorField& Sf,
        const vectorField& uS,
        const scalarField& rhoS,
        const scalarField& pS,
        label iObs,
        label iSurf,
        scalar ct
    );

    //- Calculates acoustic data used to compute acoustic pressure
    // for the specified observer
    virtual void calculateAcousticData
    (
        const vectorField& Sf,
        const vectorField& uS,
        const scalarField& rhoS,
        const scalarField& pS,
        label iObs,
        label iSurf,
        scalar ct
    ) = 0;
    
    //-
    virtual void update();
    
    //- Cleares memory from the old acoustic data. Data is considered as old
    //when it can not be used anymore by any observer because distance between
    //face and the observer is larger than c*dt, where c - sonic speed, dt -
    //time needed for an acoustic wave to get from the face to observer.
    void clearExpiredData();

    //- Dumps acoustic to the file
    void dumpAcousticData(label iObs, label iSurf);
};

}; //namespace functionObjects

}; // namespace Foam

#endif

//
//END-OF-FILE
//



/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ReadFromFile.H"
#include "FfowcsWilliamsHawkings.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ReadFromFile::ReadFromFile
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    qdsInput_(nullptr),
    fdsInput_(nullptr),
    n_blocks_to_read_(2)
{
    if (fwh_.operRegime_ == FfowcsWilliamsHawkings::DUMP_DATA)
    {
        FatalError <<
            "ReadFromFile formulation cannot work with "
            "dumpToFile regime, see the corresponding dict"
            << nl << abort(FatalError);
    }
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ReadFromFile::~ReadFromFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::ReadFromFile::calculateAcousticData
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    while(!qds_in(iObs,iSurf).eof() && !fds_in(iObs,iSurf).eof())
    {
        // Step 1. Find the maximum available time point in qds_ / fds_.
        bool t_out_of_bounds = false;
        forAll(qds_[iObs][iSurf], iFace)
        {
            if (!qds_[iObs][iSurf][iFace].first().size())
            {
                t_out_of_bounds = true;
                break;
            }
            if (qds_[iObs][iSurf][iFace].first().back() <= ct)
            {
                t_out_of_bounds = true;
                break;
            }
        }

        // Step 2. If it is larger than ct, then load data from files.
        if (t_out_of_bounds)
        {
            Info<< "Loading new acoustic data from file for time "
                << ct
                << ", qds eof = " << qds_in(iObs,iSurf).eof()
                << ", fds eof = " << fds_in(iObs,iSurf).eof()
                << endl;
            pullDataFromFiles(iObs, iSurf);
        }

        if (!t_out_of_bounds)
        {
            break;
        }
    }
}

void Foam::functionObjects::ReadFromFile::update()
{
    fwhFormulation::update();
}

void Foam::functionObjects::ReadFromFile::initialize()
{
    qdsInput_.reset(new IFstream(qdsName(), IOstreamOption::BINARY));
    fdsInput_.reset(new IFstream(fdsName(), IOstreamOption::BINARY));

    Info << "qds name = " << qdsName() << endl;
    Info << "qds opened = " << qdsInput_().opened() << endl;
    Info << "qds good = " << qdsInput_().good() << endl;
    
    //! Check for opened & good
}

Foam::IFstream&
Foam::functionObjects::ReadFromFile::qds_in(label iSurf, label iObs)
{
    return qdsInput_();
}

Foam::IFstream& Foam::functionObjects::ReadFromFile::fds_in(label iSurf, label iObs)
{
    return fdsInput_();
}

void Foam::functionObjects::ReadFromFile::pullDataFromFiles
(
    label iObs, label iSurf
)
{
    label i_block = 0;
    List<scalar> tobs_buffer (qds_[iObs][iSurf].size());
    List<scalar> xds_buffer (qds_[iObs][iSurf].size());
    // read the block
    while (n_blocks_to_read_ > i_block++)
    {
        //read t_obs, qds
        tobs_buffer = 0.0;
        xds_buffer = 0.0;
        qds_in(iObs,iSurf).readRaw(reinterpret_cast<char*>(tobs_buffer.data()), sizeof(scalar)*tobs_buffer.size());
        qds_in(iObs,iSurf).readRaw(reinterpret_cast<char*>(xds_buffer.data()), sizeof(scalar)*xds_buffer.size());
        Info << i_block << ", max/min tobs: " << max(tobs_buffer) << "/" << min(tobs_buffer) << ", sz = " << (sizeof(scalar)*tobs_buffer.size()) << endl;
        Info << i_block << ", max/min xds: " << max(xds_buffer) << "/" << min(xds_buffer) << endl;
        Info<<"tobs = " << tobs_buffer << endl;
        forAll(qds_[iObs][iSurf], iFace)
        {
            qds_[iObs][iSurf][iFace].first().append(tobs_buffer[iFace]);
            qds_[iObs][iSurf][iFace].second().append(xds_buffer[iFace]);
        }

        tobs_buffer = 0.0;
        xds_buffer = 0.0;
        fds_in(iObs,iSurf).readRaw(reinterpret_cast<char*>(tobs_buffer.data()), sizeof(scalar)*tobs_buffer.size());
        fds_in(iObs,iSurf).readRaw(reinterpret_cast<char*>(xds_buffer.data()), sizeof(scalar)*xds_buffer.size());
        forAll(fds_[iObs][iSurf], iFace)
        {
            fds_[iObs][iSurf][iFace].first().append(tobs_buffer[iFace]);
            fds_[iObs][iSurf][iFace].second().append(xds_buffer[iFace]);
        }

        // ???
        if (fds_in(iObs,iSurf).eof() || qds_in(iObs,iSurf).eof())
        {
            break;
        }
    }
}

// ************************************************************************* //

//
//END OF FILE
//
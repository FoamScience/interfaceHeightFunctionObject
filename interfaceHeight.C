#include "interfaceHeight.H"
#include "interpolation.H"
#include "midPointAndFaceSet.H"
#include "foamTime.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceHeight, 0);
    addToRunTimeSelectionTable(functionObject, interfaceHeight, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writePositions()
{
    const uniformDimensionedVectorField& g =
        mesh_.time().lookupObject<uniformDimensionedVectorField>("g");
    vector gHat = vector::zero;

    if (mag(direction_) > 0.0)
    {
        gHat = direction_/mag(direction_);
    }
    else
    {
        gHat = g.value()/mag(g.value());
    }

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    autoPtr<interpolation<scalar>>
        interpolator
        (
            interpolation<scalar>::New(interpolationScheme_, alpha)
        );

    if (Pstream::master())
    {
        //files(0) << mesh_.time().timeName() << tab;
        files(1) << mesh_.time().timeName() << tab;
    }

    forAll(locations_, li)
    {
        // Create a set along a ray projected in the direction of gravity
        const midPointAndFaceSet set
        (
            "",
            mesh_,
            ms_,
            "xyz",
            locations_[li] + gHat*mesh_.bounds().mag(),
            locations_[li] - gHat*mesh_.bounds().mag()
        );

        // Find the height of the location above the boundary
        scalar hLB = set.size() ? - gHat & (locations_[li] - set[0]) : - GREAT;
        reduce(hLB, maxOp<scalar>());

        // Calculate the integrals of length and length*alpha along the sampling
        // line. The latter is equal to the equivalent length with alpha equal
        // to one.
        scalar sumLength = 0, sumLengthAlpha = 0;
        for(label si = 0; si < set.size() - 1; ++ si)
        {
            if (set.segments()[si] != set.segments()[si+1])
            {
                continue;
            }

            const vector& p0 = set[si], p1 = set[si+1];
            const label c0 = set.cells()[si], c1 = set.cells()[si+1];
            const label f0 = set.faces()[si], f1 = set.faces()[si+1];
            const scalar a0 = interpolator->interpolate(p0, c0, f0);
            const scalar a1 = interpolator->interpolate(p1, c1, f1);

            const scalar l = - gHat & (p1 - p0);
            sumLength += l;
            sumLengthAlpha += l*(a0 + a1)/2;
        }

        reduce(sumLength, sumOp<scalar>());
        reduce(sumLengthAlpha, sumOp<scalar>());

        // Write out
        if (Pstream::master())
        {
            // Interface heights above the boundary and location
            const scalar hIB =
                liquid_ ? sumLengthAlpha : sumLength - sumLengthAlpha;
            const scalar hIL = hIB - hLB;

            // Position of the interface
            const point p = locations_[li] - gHat*hIL;

            const Foam::Omanip<int> w = setw(10u + 8u + 1u);

            //files(0) << w << hIB << w << hIL;
            files(1) << '(' << w << p.x() << w << p.y()
                << setw(10u + 8u + 0u) << p.z() << ") ";
        }
    }

    if (Pstream::master())
    {
        //files(0).endl();
        files(1).endl();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writeFileHeader(const label fid)
{
    forAll(locations_, li)
    {
        writeHeaderValue(fid, "Location " + Foam::name(li), locations_[li]);
    }

    switch (fid)
    {
        case 0:
        {
            word property = "hB";
            word value = "Interface height above the boundary";
            writeHeaderValue(fid, "hB", 
                "Interface height above the boundary");
            writeHeaderValue
            (
                fid,
                "hL",
                "Interface height above the location"
            );
            break;
        }
        case 1:
        {
            writeHeaderValue(fid, "p", "Interface position");
            break;
        }
    }

    const Foam::Omanip<int> w = setw(10u + 8u + 1u);

    forAll(locations_, li)
    {
        switch (fid)
        {
            //case 0:
            //    files(fid) << w << li << w << ' ';
            //    break;
            case 1:
                files(fid) << w << li << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(fid).endl();

    forAll(locations_, li)
    {
        switch (fid)
        {
            //case 0:
            //    files(fid) << w << "hB" << w << "hL";
            //    break;
            case 1:
                files(fid) << w << "p" << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(fid).endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceHeight::interfaceHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    mesh_(
        runTime.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", polyMesh::defaultRegion)
        )
    ),
    ms_(mesh_),
    liquid_(true),
    alphaName_("alpha"),
    interpolationScheme_("cellPoint"),
    direction_(vector::zero),
    locations_()
{
    read(dict);
    //resetNames({"height", "position"});
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceHeight::read(const dictionary& dict)
{
    dict.readIfPresent("alpha", alphaName_);
    dict.readIfPresent("liquid", liquid_);
    locations_= List<point>(dict.lookup("locations"));
    //dict.readEntry("locations", locations_);
    dict.readIfPresent("interpolationScheme", interpolationScheme_);
    dict.readIfPresent("direction", direction_);

    return true;
}


bool Foam::functionObjects::interfaceHeight::execute(const bool)
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::end()
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::write()
{
    if (Pstream::master())
    {
        const word startTimeName =
            mesh_.time().timeName(mesh_.time().startTime().value());

        wordList names = {"heights", "positions"};
        forAll(names, i)
        {
            if (!filePtr_.valid())
            {
                autoPtr<OFstream> osPtr;

                if (Pstream::master())
                {
                    const word timeName = Time::timeName(0);
                    fileName outputDir(mesh_.time().globalCaseName()/"postProcessing");
                    mkDir(outputDir);

                    word fName(names[i]);

                    // Check if file already exists
                    IFstream is(outputDir/(fName + ".dat"));
                    if (is.good())
                    {
                        fName = fName + "_" + timeName;
                    }

                    osPtr.reset(new OFstream(outputDir/(fName + ".dat")));

                    if (!osPtr->good())
                    {
                        FatalIOErrorInFunction(osPtr()) << "Cannot open file"
                            << exit(FatalIOError);
                    }

                    osPtr().setf(ios_base::scientific, ios_base::floatfield);
                    osPtr().precision(10u);
                    osPtr().width(10u + 8u);
                    filePtr_ = osPtr;
                }

                //initStream(filePtrs_[i]);
                filePtr_->setf(ios_base::scientific, ios_base::floatfield);
                filePtr_->precision(10u);
                filePtr_->width(10u + 8u);
            }
        }
    }

    writePositions();

    return true;
}


// ************************************************************************* //

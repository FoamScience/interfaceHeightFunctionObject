#ifndef interfaceHeight_H
#define interfaceHeight_H

#include "functionObject.H"
//#include "fvMeshFunctionObject.H"
//#include "logFiles.H"
#include "point.H"
#include "OFstream.H"
#include "Pstream.H"
#include "fvMesh.H"
#include "IOmanip.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class interfaceHeight Declaration
\*---------------------------------------------------------------------------*/

class interfaceHeight
:
    public functionObject
    //public logFiles
{
    // Private Data
        const fvMesh& mesh_;
        meshSearch ms_;

        //- Is the alpha field that of the liquid under the wave?
        bool liquid_;

        //- Name of the alpha field
        word alphaName_;

        //- Interpolation scheme
        word interpolationScheme_;

        //- Direction of interface motion
        vector direction_;

        //- List of locations to report the height at
        List<point> locations_;

        //- List of file pointers
        autoPtr<OFstream> filePtr_;
        //PtrList<OFstream> filePtrs_;

    // Private Member Functions

        //- Output positions
        void writePositions();


    // Private Enumerations

        //- File enumeration, not used anymore
        enum class fileID
        {
            heightFile = 0,
            positionFile = 1
        };


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label fid);
        template<class T>
        void writeHeaderValue(const label fid, const word& property, const T& value)
        {
            files(fid)  << setw(1) << '#' << setw(1) << ' '
                << setf(ios_base::left) << setw(10u + 8u - 2) << property.c_str()
                << setw(1) << ':' << setw(1) << ' ' << value << nl;
        }


public:

    //- Runtime type information
    TypeName("interfaceHeight");


    // Constructors

        //- Construct from Time and dictionary
        interfaceHeight
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        interfaceHeight(const interfaceHeight&) = delete;

        //- No copy assignment
        void operator=(const interfaceHeight&) = delete;


    //- Destructor
    virtual ~interfaceHeight() = default;


    // Member Functions

        //- Return file corresponding to enumeration
        OFstream& files(const label fid)
        {
            if (!Pstream::master())
            {
                FatalErrorInFunction
                    << "Request for files() can only be done by the master process"
                    << abort(FatalError);
            }
            //if (!filePtrs_[fid])
            //{
            //    FatalErrorInFunction
            //        << "File pointer at index " << fid << " not allocated"
            //        << abort(FatalError);
            //}

            return filePtr_();
        }


        //- Read
        virtual bool read(const dictionary&);

        //- Read
        virtual bool start(){return true;}

        //- Execute
        virtual bool execute(const bool forceWrite = false);

        //- Execute at the final time-loop
        virtual bool end();

        virtual void updateMesh(const mapPolyMesh& mpm){};
        virtual void movePoints(const pointField& points){};

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

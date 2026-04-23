#include "volFields.H"
#include "vectorField.H"
#include "vector.H"
#include "nanoflann.hpp"
#include "argList.H"
#include "IFstream.H"
#include <unordered_map>

using namespace nanoflann;
using namespace Foam;

struct FieldData
{
    word name;
    word fieldName;
    bool isVector;

    label scalarColumn;
    FixedList<label,3> vectorColumns;

    dimensionSet dims;

    std::vector<scalar> scalarValues;
    std::vector<vector> vectorValues;

    FieldData()
    :
        name(""),
        fieldName(""),
        isVector(false),
        scalarColumn(-1),
        vectorColumns({-1,-1,-1}),
        dims(dimless)
    {}
};

struct PointCloud
{
    std::vector<vector> pts;
    std::vector<FieldData> fields;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if(dim==0) return pts[idx].x();
        if(dim==1) return pts[idx].y();
        return pts[idx].z();
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree =
KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    3
>;

int main(int argc, char *argv[])
{

    argList args(argc, argv);

    #include "createTime.H"
    #include "createMesh.H"


    Info<<"Reading dictionary\n"<<endl;

    IOdictionary dict
    (
        IOobject
        (
            "mapFieldsCSVDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    fileName csvFile(dict.lookup("file"));


    dictionary interp(dict.subDict("interpolation"));

    label K = interp.lookupOrDefault<label>("neighbours", 8);
    scalar power = interp.lookupOrDefault<scalar>("power", 2.0);


    PointCloud cloud;


    // Read fields
    const dictionary& fieldsDict = dict.subDict("fields");

    forAllConstIter(dictionary, fieldsDict, iter)
    {
        FieldData field;

        field.name = iter().keyword();
        std::cout << field.name << std::endl;

        const dictionary& fDict = fieldsDict.subDict(field.name);
        Info << fDict << endl;

        word type(fDict.lookup("type"));
        field.fieldName = fDict.lookupOrDefault<word>("fieldName", field.name);

        field.dims.reset(fDict.lookup<dimensionSet>("dimensions"));
        if(type == "scalar")
        {
            field.isVector = false;
            field.scalarColumn = fDict.lookup<label>("column");
        }
        else if(type == "vector")
        {
            field.isVector = true;
            List<label> cols(fDict.lookup("columns"));
            for(label i=0;i<3;i++)
                field.vectorColumns[i] = cols[i];
        }

        cloud.fields.push_back(field);

    }


    Info<<"Reading CSV: "<<csvFile<<endl;

    std::ifstream in(csvFile);
    std::string line;
    while(std::getline(in,line))
    {
        if(line.empty() || line[0]=='#') continue;

        std::stringstream ss(line);
        std::string token;

        std::vector<std::string> tokens;

        while(std::getline(ss,token,','))
        {
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t")+1);
            tokens.push_back(token);
        }

        if(tokens.size() < 3) continue;

        scalar x = std::stod(tokens[0]);
        scalar y = std::stod(tokens[1]);
        scalar z = std::stod(tokens[2]);

        cloud.pts.emplace_back(x,y,z);


        for(auto& field : cloud.fields)
        {
            if(field.isVector)
            {
                scalar vx =
                    std::stod(tokens[field.vectorColumns[0]]);
                scalar vy =
                    std::stod(tokens[field.vectorColumns[1]]);
                scalar vz =
                    std::stod(tokens[field.vectorColumns[2]]);

                field.vectorValues.emplace_back(vx,vy,vz);
            }
            else
            {
                scalar v =
                    std::stod(tokens[field.scalarColumn]);

                field.scalarValues.push_back(v);
            }
        }
    }


    Info<<"Building KDTree"<<endl;

    KDTree tree(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();


    List<autoPtr<volScalarField>> scalarFields;
    List<autoPtr<volVectorField>> vectorFields;


    for(auto& f : cloud.fields)
    {
        if(f.isVector)
        {
            vectorFields.append
            (
                autoPtr<volVectorField>
                (
                    new volVectorField
                    (
                        IOobject
                        (
                            f.fieldName,
                            runTime.time().constant(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedVector
                        (
                            "zero",
                            f.dims,
                            vector(0,0,0)
                        )
                    )
                )
            );
        }
        else
        {
            std::cout << "Appending scalar fields " << std::endl;
            scalarFields.append
            (
                autoPtr<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            f.fieldName,
                            runTime.time().constant(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "zero",
                            f.dims,
                            0
                        )
                    )
                )
            );
        }
    }


    std::vector<size_t> idx(K);
    std::vector<double> dist2(K);


    forAll(mesh.C(), cellI)
{
    vector pt = mesh.C()[cellI];

    double query[3] = {pt.x(),pt.y(),pt.z()};

    nanoflann::KNNResultSet<double> resultSet(K);
    resultSet.init(idx.data(),dist2.data());

    tree.findNeighbors
    (
        resultSet,
        query,
        nanoflann::SearchParameters(10)
    );

    label sId = 0;
    label vId = 0;

    for(auto& f : cloud.fields)
    {
        scalar sumW = 0;

        vector sumV(0,0,0);
        scalar sumS = 0;

        // for deviation check
        scalar minVal = GREAT;
        scalar maxVal = -GREAT;

        for(size_t i=0;i<resultSet.size();i++)
        {
            scalar d = std::sqrt(dist2[i]) + SMALL;
            scalar w = 1.0/std::pow(d,power);

            sumW += w;

            if(f.isVector)
            {
                const vector& v = f.vectorValues[idx[i]];
                sumV += w * v;

                scalar magV = mag(v);
                minVal = std::min(minVal, magV);
                maxVal = std::max(maxVal, magV);
            }
            else
            {
                scalar v = f.scalarValues[idx[i]];
                sumS += w * v;

                minVal = std::min(minVal, v);
                maxVal = std::max(maxVal, v);
            }
        }

        scalar deviation = 0.0;
        if(f.isVector)
        {
            vector interp = sumV / (sumW + SMALL);
            scalar meanMag = mag(interp);

            if(meanMag > SMALL)
                deviation = (maxVal - minVal) / meanMag;

            if(sumW > SMALL && deviation < 0.4)
            {
                (*vectorFields[vId])[cellI] = interp;
            }
            else
            {
                // fallback = nearest neighbor
                (*vectorFields[vId])[cellI] =
                    f.vectorValues[idx[0]];
            }

            vId++;
        }
        else
        {
            scalar interp = sumS / (sumW + SMALL);

            if(std::abs(interp) > SMALL)
                deviation = (maxVal - minVal) / std::abs(interp);

            if(sumW > SMALL && deviation < 0.4)
            {
                (*scalarFields[sId])[cellI] = interp;
            }
            else
            {
                (*scalarFields[sId])[cellI] =
                    f.scalarValues[idx[0]];
            }

            sId++;
        }
    }
}


    Info<<"Writing fields"<<endl;

    for(auto& f : scalarFields)
        f().write();

    for(auto& f : vectorFields)
        f().write();


    Info<<"Done"<<endl;

    return 0;
}
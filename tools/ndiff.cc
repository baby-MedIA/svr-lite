// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR

#define SVR_TOOL
#include "svr/Common.h"
#include "svr/Profiling.h"

using namespace std;
using namespace svr;
using namespace boost::program_options;

// options_description getOptions() {
//     options_description opts("Options");
//     return opts;
// }

// void printUsage() {
// }

template<typename T>
inline bool CompareAttribute(const T& a, const T& b) {
    if(a==b)
        return true;
    if (std::isnan(a) && std::isnan(b))
        return true;
    return false;
}

template bool CompareAttribute<int>(const int& a, const int& b);
template bool CompareAttribute<float>(const float& a, const float& b);
template bool CompareAttribute<double>(const double& a, const double& b);


bool EqualImages(const RealImage& imgA, const RealImage& imgB, float dataTolerance, bool verbose, bool createMap) {

    bool match = true;
    const ImageAttributes& imgAattr = imgA.Attributes();
    const ImageAttributes& imgBattr = imgB.Attributes();
    RealImage map(imgAattr);

    // Compare attributes - these must match
    if (!CompareAttribute(imgAattr._x,imgBattr._x))   {
        cout << "_x differs: " << imgAattr._x << " vs " << imgBattr._x << endl;
        match = false;
    }
    if (!CompareAttribute(imgAattr._y , imgBattr._y))   {
        cout << "_y differs: " << imgAattr._y << " vs " << imgBattr._y << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._z , imgBattr._z))   {
        cout << "_z differs: " << imgAattr._z << " vs " << imgBattr._z << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._t , imgBattr._t))   {
        cout << "_t differs: " << imgAattr._t << " vs " << imgBattr._t << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._dx , imgBattr._dx))   {
        cout << "_dx differs: " << imgAattr._dx << " vs " << imgBattr._dx << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._dy , imgBattr._dy))   {
        cout << "_dy differs: " << imgAattr._dy << " vs " << imgBattr._dy << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._dz , imgBattr._dz))   {
        cout << "_dz differs: " << imgAattr._dz << " vs " << imgBattr._dz << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._dt , imgBattr._dt))   {
        cout << "_dt differs: " << imgAattr._dt << " vs " << imgBattr._dt << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._xorigin , imgBattr._xorigin))   {
        cout << "_xorigin differs: " << imgAattr._xorigin << " vs " << imgBattr._xorigin << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._yorigin , imgBattr._yorigin))   {
        cout << "_yorigin differs: " << imgAattr._yorigin << " vs " << imgBattr._yorigin << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._zorigin , imgBattr._zorigin))   {
        cout << "_zorigin differs: " << imgAattr._zorigin << " vs " << imgBattr._zorigin << endl;
        match =  false;
    }
    if (!CompareAttribute(imgAattr._torigin , imgBattr._torigin))   {
        cout << "_torigin differs: " << imgAattr._torigin << " vs " << imgBattr._torigin << endl;
        match =  false;
    }
    
    for (int i = 0; i < 3; i++){
        if (match && imgAattr._xaxis[i] != imgBattr._xaxis[i])   {
            cout << "_xaxis differs: " << i << "->" << imgAattr._xaxis[i] << " vs " << imgBattr._xaxis[i] << endl;
            match =  false;
        }
        if (match && imgAattr._yaxis[i]  != imgBattr._yaxis[i] )   {
            cout << "_yaxis differs: " << i << "->" << imgAattr._yaxis[i] << " vs " << imgBattr._yaxis[i] << endl;
            match =  false;
        }
        if (match && imgAattr._zaxis[i]  != imgBattr._zaxis[i] )   {
            cout << "_zaxis differs: " << i << "->" << imgAattr._zaxis[i] << " vs " << imgBattr._zaxis[i] << endl;
            match =  false;
        }
    }

    if (match && !CompareAttribute(imgA.NumberOfVoxels() , imgB.NumberOfVoxels()))   {
        cout << "Number of Voxels differs: " << imgA.NumberOfVoxels() 
                << " vs " << imgB.NumberOfVoxels() << endl;
        match =  false;
    }

    // compare voxel intensities
    int count = 0;
    for (int t = 0; t < imgAattr._t; ++t)
            for (int z = 0; z < imgAattr._z; ++z)
                for (int y = 0; y < imgAattr._y; ++y)
                    for (int x = 0; x < imgAattr._x; ++x) {
                        float aVal = imgA.Get(x,y,z,t);
                        float bVal = imgB.Get(x,y,z,t);
                        float diff = aVal-bVal;
                        if(abs(diff) > dataTolerance) {
                            cout << setprecision(9);
                            if (verbose || match){
                                cout << "Image voxel difference at x=" << x << " y=" << y << " z=" << z << " t=" << t;
                                cout << ": Image A intensity=" << aVal << " vs Image B intensity=" << bVal << endl;
                            }
                            match =  false;
                            ++count;
                            if(createMap) map(x,y,z,t)=1;
                        }
                        else {
                            if(createMap) map(x,y,z,t)=0;
                        }
                    }
    if (count > 0) cout << "Number of different voxels = " << count << endl;
    if(createMap) map.Write("map.nii.gz");
    return match;
}

int main(int argc, char **argv) {

    // Define option parameters
    options_description opts("Options");
    variables_map vm;

    options_description reqOpts;
    string target;
    string source;
    double tolerance;
    bool verbose;
    bool map;
    reqOpts.add_options()
        ("img1", value<string>(&target)->required(), "Image 1")
        ("img2", value<string>(&source)->required(), "Image 2");
        
    options_description allOpts("Allowed options");
    allOpts.add(reqOpts).add(opts).add_options() 
        ("tol", value<double>(&tolerance), "Difference tolerance [Default: 0 (no tolerance)]")
        ("map", bool_switch(&map), "Output a binary diff map to specified path")
        ("verbose", bool_switch(&verbose), "Verbose output");
    
    positional_options_description posOpts;
    posOpts.add("img1", 1).add("img2", 1);

    // Read input parameters
    try {
        store(command_line_parser(argc, argv).options(allOpts).positional(posOpts)
            .style(command_line_style::unix_style | command_line_style::allow_long_disguise).run(),vm);
        notify(vm);
    } catch (error& e) {
        // Delete -- from the argument name in the error message
        string err = e.what();
        cerr << "Argument parsing error: " << err << "\n\n";
        //printUsage();
        return 1;
    }

    InitializeIOLibrary();

    UniquePtr<ImageReader> image_reader;
    UniquePtr<BaseImage> tmp_image;

    // Read image1
    image_reader.reset(ImageReader::TryNew(target.c_str()));
    tmp_image.reset(image_reader->Run());
    RealImage targetImg;
    targetImg = *tmp_image;

    // Read img2
    image_reader.reset(ImageReader::TryNew(source.c_str()));
    tmp_image.reset(image_reader->Run());
    RealImage sourceImg;
    sourceImg = *tmp_image;

    //Compare target and simulated target within mask region
    bool result = EqualImages(targetImg, sourceImg, tolerance, verbose, map);
    cout << "Match=" << result << endl;
}
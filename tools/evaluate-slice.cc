// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR

#define SVR_TOOL
#include "svr/Profiling.h"
#include "reconstruct.h"

using namespace std;
using namespace svr;
using namespace boost::program_options;

options_description getOptions() {
    options_description opts("Options");
    return opts;
}

void printUsage() {
}

// transform mask to the specified source image
void TransformImage(const RealImage& target, RealImage& source, const RigidTransformation& transformation, int sourcePadding) {
    //transform mask to the space of image
    unique_ptr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(Interpolation_NN));
    ImageTransformation imagetransformation;
    RealImage m(target.Attributes());

    imagetransformation.Input(&source);
    imagetransformation.Transformation(&transformation);
    imagetransformation.Output(&m);
    //target contains zeros and ones image, need padding -1
    imagetransformation.TargetPaddingValue(-1);
    //need to fill voxels in target where there is no info from source with zeroes
    imagetransformation.SourcePaddingValue(sourcePadding);
    imagetransformation.Interpolator(interpolator.get());
    imagetransformation.Run();
    source = move(m);
}

/**
 * Calculate Normalised Cross Coefficient
 * as:
 * sum((targetVoxel - mean_targetVoxel) * (sourceVoxel - mean_sourceVoxel)) /
 * sqrt( sum((targetVoxel - mean_targetVoxel)^2) * sum((sourceVoxel - mean_sourceVoxel)^2)) )
 */
double CalculateNCC(const RealImage& target, RealImage& source, RealImage& mask){
        double sumTargetVal=0, countTargetVal=0, sumSourceVal=0, countSourceVal=0;
        double targetMean = 0, sourceMean = 0;

        double diff_sum = 0;
        double targetNormSq = 0;
        double sourceNormSq = 0;

        // Calculate means:
        #pragma omp parallel for reduction(+: sumTargetVal, countTargetVal, sumSourceVal, countSourceVal) collapse(3)
        for (int l=0; l<target.GetT(); l++) {
            for (int j=0; j<target.GetY(); j++) {
                for (int i=0; i<target.GetX(); i++) {
                    if(mask(i,j,0,l) > 0) {
                        double targetValue = target(i,j,0,l);
                        double sourceValue = source(i,j,0,l);

                        sumTargetVal += targetValue;
                        countTargetVal++;
                        sumSourceVal += sourceValue;
                        countSourceVal++;
                    }
                }
            }
        }
        targetMean = sumTargetVal / countTargetVal;
        sourceMean = sumSourceVal / countSourceVal;

        // Calculate NCC:
        #pragma omp parallel for reduction(+: diff_sum, targetNormSq, sourceNormSq) collapse(3)
        for (int l=0; l<target.GetT(); l++) {
            for (int j=0; j<target.GetY(); j++) {
                for (int i=0; i<target.GetX(); i++) {
                    if(mask(i,j,0,l) > 0) {
                        double targetValue = target(i,j,0,l);
                        double sourceValue = source(i,j,0,l);

                        diff_sum += (targetValue - targetMean) * (sourceValue - sourceMean);
                        targetNormSq += pow(targetValue - targetMean, 2);
                        sourceNormSq += pow(sourceValue - sourceMean, 2);
                    }
                }
            }
        }

        double ncc = 0;
        if (targetNormSq * sourceNormSq > 0)
            ncc = diff_sum / sqrt(targetNormSq * sourceNormSq);
        
        return ncc;
}


int main(int argc, char **argv) {

    // Define option parameters
    options_description opts = getOptions();
    variables_map vm;

    options_description reqOpts;
    string target;
    string source;
    string dof;
    string mask;
    reqOpts.add_options()
        ("target", value<string>(&target)->required(), "Target 2D image slice")
        ("source", value<string>(&source)->required(), "Source 3D image volume")
        ("dof", value<string>(&dof)->required(), "Transformation to map 2D slice to 3D volume.")
        ("mask", value<string>(&mask)->required(), "Mask region of interest in 3D volume.");
        
    options_description allOpts("Allowed options");
    allOpts.add(reqOpts).add(opts);
    
    positional_options_description posOpts;
    posOpts.add("target", 1).add("source", 1).add("dof", 1).add("mask", 1);

    // Read input parameters
    try {
        store(command_line_parser(argc, argv).options(allOpts).positional(posOpts)
            .style(command_line_style::unix_style | command_line_style::allow_long_disguise).run(),vm);
        notify(vm);
    } catch (error& e) {
        // Delete -- from the argument name in the error message
        string err = e.what();
        cerr << "Argument parsing error: " << err << "\n\n";
        printUsage();
        return 1;
    }

    InitializeIOLibrary();

    UniquePtr<ImageReader> image_reader;
    UniquePtr<BaseImage> tmp_image;

    // Read target (slice) image
    image_reader.reset(ImageReader::TryNew(target.c_str()));
    tmp_image.reset(image_reader->Run());
    RealImage targetImg;
    targetImg = *tmp_image;

    // Read source (volume) image
    image_reader.reset(ImageReader::TryNew(source.c_str()));
    tmp_image.reset(image_reader->Run());
    RealImage sourceImg;
    sourceImg = *tmp_image;

    // Read transformation (dof) of target->source
    Transformation *transformation = Transformation::New(dof.c_str());
    unique_ptr<RigidTransformation> rigidTransf(dynamic_cast<RigidTransformation*>(transformation));

    // Read mask image
    image_reader.reset(ImageReader::TryNew(mask.c_str()));
    tmp_image.reset(image_reader->Run());
    RealImage maskImg;
    maskImg = *tmp_image;

    //Create simulated slice images from source volume and mask volume
    //rigidTransf->Invert();
    TransformImage(targetImg, sourceImg, *rigidTransf, -1);
    TransformImage(targetImg, maskImg, *rigidTransf, 0);

    sourceImg.Write("debug_source.nii.gz");
    targetImg.Write("debug_target.nii.gz");
    maskImg.Write("debug_mask.nii.gz");

    //Compare target and simulated target within mask region
    double ncc = CalculateNCC(targetImg, sourceImg, maskImg);
    cout << "NCC=" << ncc << endl;
}
// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR
#include "svr/Evaluation.h"
#include "svr/Profiling.h"

using namespace std;

namespace svr {

    Evaluation::Evaluation(const bool debug, const bool verbose, const bool profile) 
    {
        _verbose = verbose;
        _debug = debug || verbose;
        _profile = profile;

        ConfigureLogging("evaluation");
    }

    int Evaluation::run(const string& stackFile,
                        double& thickness,
                        const string& maskFile,
                        const string& reconVolume,
                        const bool stackIntensityMatching,
                        const int& iterations) {

        START_TIMING();

        // Read Stacks
        int tSize;
        RealImage rawStack;
        ReadStack(stackFile, rawStack, tSize);

        // Read Masks
        RealImage mask;
        ReadMask(maskFile, mask);

        // Calculate slice thickness
        CalculateSliceThickness(rawStack, thickness); 

        // Read volume
        RealImage volume;
        int tVolSize;
        ReadVolume(reconVolume, volume, tVolSize);
        if (_debug) volume.Write("volume.nii.gz");

        // Extract mask from volume padding
        RealImage volumeMask;
        ExtractVolumeMask(volume, volumeMask);
        if (_debug) volumeMask.Write("volume-mask.nii.gz");

        // Create reconstruction model
        ReconstructionModel reconstruction;
        reconstruction.SetReconstructed(volume);
        reconstruction.SetMask(volumeMask,0);
        if (_profile) reconstruction.ProfileOn();
        else reconstruction.ProfileOff();

        // Register stack to volume
        vector<RigidTransformation> transformations(1);
        vector<RealImage> stacks  = { rawStack };
        reconstruction.StackRegistrationsWithTemplate({stacks}, transformations, volume);
        ReconstructionModel::SaveStackTransformation(transformations[0], 0);


        // Stack intensity matching
        ReconstructionParameters parameters; // use only for default values
        if (stackIntensityMatching) {
            reconstruction.MatchStackIntensitiesWithMasking(stacks, transformations, parameters.averageValue);  
            if (_debug) reconstruction.OutputStackFactor();
        }

        // Evaluation stack list
        // Create slices
        reconstruction.CreateSlicesAndTransformations(stacks, transformations, {thickness});
        reconstruction.CreateSliceMasks(); // uses volume mask

        reconstruction.InitializeEM(); //Initialise _parameters structures for EM

        // Register slices to volume
        reconstruction.SliceToVolumeRegistration();
        reconstruction.SaveTransformations();

        //Initialise values of weights, scales and bias fields
        reconstruction.InitializeEMValues();       

        // Create coefficient matrix
        reconstruction.CoeffInit(); // uses volume mask

        // Simulate slices
        // Input number of iterations, output similarity measure between aquired and simulated.
        reconstruction.SimulateSlices(); // uses volume mask. Needs intensity matching. Output bias and scale and monitor changes to determine iterations.
        if(_debug)reconstruction.SaveSimulatedSlices();

        // Simulate slices
        // Input number of iterations, output similarity measure between aquired and simulated.
        reconstruction.SimulateSlices(); // uses volume mask. Needs intensity matching. Output bias and scale and monitor changes to determine iterations.
        if (_debug) reconstruction.SaveSimulatedSlices();
        reconstruction.SaveSimulatedSlicesAsStack(stacks,0);

        // Intensity matching iterations.
        logger::debug("iterations=" + to_string(iterations));
        bool converged = false;
        int count = 0;
        vector<double> consistency_history;
        int limit = 2000;
        int range = 20;
        double threshold = 0.1;
        while(!converged) {
            ++count;

            // Apply Bias and Scale
            reconstruction.Bias();   
            reconstruction.Scale();
            
            reconstruction.SaveBiasFieldsAsStack(stacks, count);

            // Calculate consistency
            double consistency = reconstruction.Consistency(true);
            if (_debug) {
                std::cout << std::fixed;
                std::cout << std::setprecision(3);
                cout << count <<  ": consistency=" << consistency << endl;
            }
            consistency_history.push_back(consistency);

            // Detect convergence
            if (iterations >= 0) { // observe user-defined iterations if defined
                if (count >= iterations) {
                    converged = true;
                }
            }
            else if (consistency_history.size() > 2) { // otherwise auto-detect convergence
                int begin = 1;
                int end = consistency_history.size()-1;
                if (consistency_history.size() > range)
                    begin = consistency_history.size() - range;

                double sum_diffs;
                double count_diffs;
                for (int idx = begin; idx < end; ++idx){
                    ++count_diffs;
                    double c1  = consistency_history[idx-1];
                    double c2 = consistency_history[idx];
                    sum_diffs += abs(c1 - c2);
                }
                double avg = sum_diffs/count_diffs;

                if (avg < threshold || count > limit) {
                    cout << "Converged count="<< count<< " average_diff=" << avg << endl;
                    converged = true;
                }
            }
        }

        // Simulate slices
        if (_debug) reconstruction.SaveSimulatedSlices();
        reconstruction.SaveSimulatedSlicesAsStack(stacks,1);

        // Get simulated stacks
        vector<RealImage> simStacks;
        reconstruction.GetSimulatedSlicesAsStacks(stacks, simStacks);
        reconstruction.SaveSimulatedSlicesAsStack(stacks, 2);
        reconstruction.SaveScaledSimulatedSlicesAsStack(stacks);
        vector<RealImage> sliceBias;
        reconstruction.GetBiasFieldsAsStack(stacks, sliceBias);
        sliceBias[0].Write("slice0_bias.nii.gz");
        reconstruction.SaveBiasFieldsAsStack(stacks);
        reconstruction.SaveAsStack(sliceBias,stacks, "finalBias");
        vector<double> sliceScale = reconstruction.GetSliceScale();

        double recon_consistency = reconstruction.Consistency(true);
        double recon_nrmse = reconstruction.CalculateNRMSE(true);
        double recon_ncc = reconstruction.CalculateNCC(true);

        // Calculate similarity of stack to simulated stack
        rawStack.Write("stack.nii.gz");
        mask.Write("cardiacMask.nii.gz");
        transformations[0].Write("transformation.nii.gz");

        // Calculate similarity measures
        stringstream ss;
        ss << "Stack, NCC, NRMSE" << endl;
        ss << rawStack << ",";
        double ncc = CalculateNCC(stacks[0], simStacks[0], mask, sliceBias[0], sliceScale);
        ss << ncc << ",";
        double nrmse = CalculateNRMSE(stacks[0], simStacks[0], mask, sliceBias[0], sliceScale);
        ss << nrmse << endl;

        logger::info(ss.str());

        END_TIMING("Run Evaluation");
        return 0;
    }

    /** 
     * Read Stacks
     */
    void Evaluation::ReadStack(const string &stackFile,
                                RealImage &stack,
                                int &tSize)
    {
        START_TIMING();
        tSize = 0;
        UniquePtr<ImageReader> image_reader;
        UniquePtr<BaseImage> tmp_image;
        logger::debug("Reading stack " + stackFile);
        image_reader.reset(ImageReader::TryNew(stackFile.c_str()));
        tmp_image.reset(image_reader->Run());
        stack = *tmp_image;
        string filename = "stack_input.nii.gz";
        if (_debug) cout << "Write " << filename << endl;
        if (_debug) stack.Write(filename.c_str());
        if (stack.GetT() > tSize) 
            tSize = stack.GetT();
        logger::debug("Max stack time points: " + to_string(tSize));
        END_TIMING("ReadStacks");
    }

    /**
     * Read stack masks
     */
    void Evaluation::ReadMask(const string &maskFile, RealImage &mask)
    {
        START_TIMING();
        UniquePtr<ImageReader> image_reader;
        UniquePtr<BaseImage> tmp_image;
        logger::debug("Reading stack " + maskFile);
        RealImage stack;
        image_reader.reset(ImageReader::TryNew(maskFile.c_str()));
        tmp_image.reset(image_reader->Run());
        mask = *tmp_image;
        if (_debug) stack.Write("stack_mask.nii.gz");
        END_TIMING("ReadMasks");
    }

    /**
     * Ensure the slice thickness vector is populated using the stack z-dimension resolution if not provided.
     */
    void Evaluation::CalculateSliceThickness(RealImage &stack, double &thickness)
    {
        START_TIMING();
        stringstream ss;
        ss << "Slice thickness is:";
        if (thickness < 0)
        {
            // Thickness hasn't been provided
            double dx, dy, dz;
            stack.GetPixelSize(&dx, &dy, &dz);
            thickness = dz;
            ss << " " << thickness;
        }
        else
        {
            // Thickness is provided by user
            ss << thickness << " ";
        }
        logger::debug(ss.str());
        END_TIMING("CalculateSliceThickness");
    }

    void Evaluation::ReadVolume(const string &volumeFile,
                                RealImage &volume,
                                int &tSize)
    {
        START_TIMING();
        tSize = 0;
        UniquePtr<ImageReader> image_reader;
        UniquePtr<BaseImage> tmp_image;
 
        logger::debug("Reading volume " + volume);
        image_reader.reset(ImageReader::TryNew(volumeFile.c_str()));
        tmp_image.reset(image_reader->Run());
        volume = *tmp_image;
        tSize = volume.GetT();
        
        logger::debug("Max stack time points: " + to_string(tSize));
        END_TIMING("ReadVolume");
    }

    /**
     * Recreate volume mask from volume image using the volume padding value.
     */
    void Evaluation::ExtractVolumeMask(const RealImage& volume, RealImage& mask){
        START_TIMING();
        ImageAttributes attr = volume.Attributes();
        attr._t = 1;
        mask.Initialize(attr);

        // Calculate padding from min value (either -1 or 0).
        RealPixel min, max;
        volume.GetMinMax(min, max);
        int padding = 0;
        if (min < -0.1)
            padding = -1;

        #pragma omp parallel for
        for (int z = 0; z < volume.GetZ(); ++z) {
            for (int y = 0; y < volume.GetY(); ++y) {
                for (int x = 0; x < volume.GetX(); ++x) {
                    if (volume(x,y,z,0) == padding)
                        mask(x,y,z) = 0;
                    else 
                        mask(x,y,z) = 1;
                }
            }
        }
        END_TIMING("ExtractVolumeMask");
    }

    /**
     * Calculate stack Normalised Cross Coefficient
     * as:
     * sum((targetVoxel - mean_targetVoxel) * (sourceVoxel - mean_sourceVoxel)) /
     * sqrt( sum((targetVoxel - mean_targetVoxel)^2) * sum((sourceVoxel - mean_sourceVoxel)^2)) )
     */
    double Evaluation::CalculateNCC(const RealImage& original, 
                                    const RealImage& simulated, 
                                    const RealImage& masks,
                                    const RealImage& bias,
                                    const vector<double>& scale){
        START_TIMING();

        double sumTargetVal=0, countTargetVal=0, sumSourceVal=0, countSourceVal=0;
        double targetMean = 0, sourceMean = 0;

        double diff_sum = 0;
        double targetNormSq = 0;
        double sourceNormSq = 0;

        // If stack masks provided then use them, otherwise use volume padding as mask.
        bool useMask = true;

        // Calculate means:
        #pragma omp parallel for reduction(+: sumTargetVal, countTargetVal, sumSourceVal, countSourceVal)
        for (int l=0; l<original.GetT(); l++) {
            for (int sliceNum=0; sliceNum<original.GetZ(); sliceNum++) {
                for (int j=0; j<original.GetY(); j++) {
                    for (int i=0; i<original.GetX(); i++) {
                        if( (useMask && masks(i,j,sliceNum) > 0) || (!useMask && simulated(i,j,sliceNum) > 0)) {
                            double targetValue = original(i,j,sliceNum,l);
                            double sourceValue = simulated(i,j,sliceNum,l) * exp(bias(i, j, sliceNum)) / scale[sliceNum];

                            sumTargetVal += targetValue;
                            countTargetVal++;
                            sumSourceVal += sourceValue;
                            countSourceVal++;
                        }
                    }
                }
            }
        }
        targetMean = sumTargetVal / countTargetVal;
        sourceMean = sumSourceVal / countSourceVal;

        // Calculate NCC:
        #pragma omp parallel for reduction(+: diff_sum, targetNormSq, sourceNormSq)
        for (int l=0; l<original.GetT(); l++) {
            for (int sliceNum=0; sliceNum<original.GetZ(); sliceNum++) {
                for (int j=0; j<original.GetY(); j++) {
                    for (int i=0; i<original.GetX(); i++) {
                        if( (useMask && masks(i,j,sliceNum) > 0 && simulated(i,j,sliceNum) > 0)  // if using stack masks: must be inside stack mask & padding
                        || (!useMask && simulated(i,j,sliceNum) > 0)) {                        // if using just volume padding: must be inside padding
                            double origValue = original(i,j,sliceNum,l);
                            double simValue = simulated(i,j,sliceNum,l) * exp(bias(i, j, sliceNum)) / scale[sliceNum];

                            diff_sum += (origValue - targetMean) * (simValue - sourceMean);
                            targetNormSq += pow(origValue - targetMean, 2);
                            sourceNormSq += pow(simValue - sourceMean, 2);
                        }
                    }
                }
            }
        }

        double ncc = 0;
        if (targetNormSq * sourceNormSq > 0)
            ncc = diff_sum / sqrt(targetNormSq * sourceNormSq);
            
        END_TIMING("CalculateNCC");
        return ncc;
    }

    /**
     * Calculate stack NRMSE
     * as
     * squareRoot ( sum(targetVoxel-sourceVoxel) / numberOfTargetVoxels ) ) / targetMean
     */
    double Evaluation::CalculateNRMSE(const RealImage& original, 
                                      const RealImage& simulated, 
                                      const RealImage& masks,
                                      const RealImage& bias,
                                      const vector<double>& scale) {
        START_TIMING();
        bool useMask = true;

        double sumDiff = 0, sumValue = 0, num = 0;

        // voxels
        #pragma omp parallel for reduction(+: sumDiff, sumValue, num)
        for (int l=0; l<original.GetT(); l++) {
            for (int sliceNum=0; sliceNum<original.GetZ(); sliceNum++) {
                for (int j=0; j<original.GetY(); j++) {
                    for (int i=0; i<original.GetX(); i++) {
                        if( (useMask && masks(i,j,sliceNum) > 0 && simulated(i,j,sliceNum) > 0)  // if using stack masks: must be inside stack mask & padding
                        || (!useMask && simulated(i,j,sliceNum) > 0)) {                        // if using just volume padding: must be inside padding
                            double origValue = original(i,j,sliceNum,l);
                            double simValue = simulated(i,j,sliceNum,l) * exp(bias(i, j, sliceNum)) / scale[sliceNum];
                            double diff = origValue - simValue; 
                            sumValue += origValue;
                            sumDiff += diff*diff;
                            num++;
                        }
                    }
                }
            }
        }

        double mean = sumValue/num;
        double nrmse = sqrt(sumDiff/num) / mean;

        END_TIMING("CalculateNRMSE");
        return nrmse;
    }
}

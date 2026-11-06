// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR
#include "svr/ReconstructionModel.h"
#include "svr/Profiling.h"
#include "svr/Parallel.h"
#include "svr/json.hpp"

using json = nlohmann::json;
using std::vector;

namespace svr {

    ReconstructionModel::ReconstructionModel() {
        _number_of_slices_org = 0;
        _average_thickness_org = 0;
        _step = 0.0001;
        _debug = false;
        _verbose = false;
        _quality_factor = 2;
        _sigma_bias = 12;
        _sigma_s = 0.025;
        _sigma_s2 = 0.025;
        _mix_s = 0.9;
        _mix = 0.9;
        _delta = 1;
        _lambda = 0.1;
        _alpha = (0.05 / _lambda) * _delta * _delta;
        _low_intensity_cutoff = 0.01;
        _nmi_bins = -1;
        _local_SSIM_threshold = 0.4;
        _local_SSIM_window_size = 20;
        _template_created = false;
        _have_mask = false;
        _global_bias_correction = false;
        _adaptive = false;
        _robust_slices_only = false;
        _recon_type = _3D;
        _template_flag = false;
        _no_sr = false;
        _masked_stacks = false;
        _no_masking_background = false;
        _no_ts = false;
        _global_NCC_threshold = 0.5;
        _current_iteration = 0;
    }

    // Pure virtual destructor for abstract base class
    ReconstructionModel::~ReconstructionModel() {}

//     //-------------------------------------------------------------------

    // Create average of all stacks based on the input transformations
    RealImage ReconstructionModel::CreateAverage(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations) {
        if (!_template_created)
            throw runtime_error("Please create the template before calculating the average of the stacks.");

        START_TIMING();
        const RealImage average = Utility::CreateAverage(this, stacks, stack_transformations);
        END_TIMING("CreateAverage");

        return average;
    }

    // Set given mask to the reconstruction object
    void ReconstructionModel::SetMask(RealImage &mask, double sigma, double threshold) {
        if (!_template_created)
            throw runtime_error("Please create the template before setting the mask, so that the mask can be resampled to the correct dimensions.");

        ImageAttributes attr = _reconstructed.Attributes();
        attr._t = 1; // Mask is 3D.
        _reconstructed_mask.Initialize(attr);

        //if sigma is nonzero first smooth the mask
        if (sigma > 0) {
            //blur mask
            GaussianBlurring<RealPixel> gb(sigma);

            gb.Input(&mask);
            gb.Output(&mask);
            gb.Run();

            //binarise mask
            mask = CreateMask(mask, threshold);
        }

        //resample the mask according to the template volume using identity transformation
        RigidTransformation transformation;
        ImageTransformation imagetransformation;

        //GenericNearestNeighborExtrapolateImageFunction<RealImage> interpolator;
        InterpolationMode interpolation = Interpolation_NN;
        unique_ptr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(interpolation));

        imagetransformation.Input(&mask);
        imagetransformation.Transformation(&transformation);
        imagetransformation.Output(&_reconstructed_mask);
        //target is zero image, need padding -1
        imagetransformation.TargetPaddingValue(-1);
        //need to fill voxels in target where there is no info from source with zeroes
        imagetransformation.SourcePaddingValue(0);
        imagetransformation.Interpolator(interpolator.get());
        imagetransformation.Run();

        //set flag that mask was created
        _have_mask = true;
        _evaluation_mask = _reconstructed_mask;

        // compute mask volume
        double vol = 0;
        RealPixel *pm = _evaluation_mask.Data();
        #pragma omp parallel for reduction(+: vol)
        for (int i = 0; i < _reconstructed_mask.NumberOfVoxels(); i++)
            if (pm[i] > 0.1)
                vol++;

        vol *= _reconstructed.GetXSize() * _reconstructed.GetYSize() * _reconstructed.GetZSize() / 1000;

        logger::debug((boost::format("ROI volume :  %1% cc") % vol).str());

        if (_debug)
            SaveReconImageFile(_reconstructed_mask,"mask");
    }

    //-------------------------------------------------------------------

    // // Get 3D version of 4D stacks
    void ReconstructionModel::Get3DStacks(const Array<RealImage>& stacks, Array<RealImage>& sim3dStacks) {
        sim3dStacks.reserve(stacks.size());
        for (int i = 0; i < stacks.size(); ++i) {
            int x = stacks[i].Attributes()._x;
            int y = stacks[i].Attributes()._y;
            int z = stacks[i].Attributes()._z;
            int t = 0;
            sim3dStacks.push_back(stacks[i].GetRegion(0, 0, 0, t, x, y, z, t+1));
        }
    }

    void ReconstructionModel::StackRegistrations(const Array<RealImage>& stacks, 
                                                Array<RigidTransformation>& stack_transformations){
            RealImage target = _reconstructed;
            StackRegistrations3D(stacks, stack_transformations, target);
    }
    void ReconstructionModel::StackRegistrations(Array<RealImage>& stacks, 
                                                Array<RigidTransformation>& stack_transformations, 
                                                const int templateNumber){
        // If 4D convert stacks to 3D
        if(stacks[templateNumber].GetT() > 1) {
            logger::debug("Converting 4D stacks to 3D stacks for registration");
            Array<RealImage> sim3dStacks;
            Get3DStacks(stacks, sim3dStacks);
            RealImage target = sim3dStacks[templateNumber];
            StackRegistrations3D(sim3dStacks, stack_transformations, target);
        }
        else {
            RealImage target = stacks[templateNumber];
            StackRegistrations3D(stacks, stack_transformations, target);
        }

    }
    void ReconstructionModel::StackRegistrationsWithTemplate(const Array<RealImage>& stacks, 
                                                Array<RigidTransformation>& stack_transformations, 
                                                RealImage input_template){
                                     
        RealImage inputTemplate3D;

        // If 4D convert template to 3D
        if (input_template.GetT() > 1) {
            inputTemplate3D = input_template.GetRegion(0,0,0,0, input_template.GetX(), input_template.GetY(), input_template.GetZ(), 1);
        }
        else {
            inputTemplate3D = input_template;
        }

        // if 4D convert stacks to 3D
        if(stacks.size() > 0 && stacks[0].GetT() > 1) {
            logger::debug("Converting 4D stacks to 3D stacks for registration");
            Array<RealImage> sim3dStacks;
            Get3DStacks(stacks, sim3dStacks);
            StackRegistrations3D(sim3dStacks, stack_transformations, inputTemplate3D);
        }
        else {
            StackRegistrations3D(stacks, stack_transformations, inputTemplate3D);
        }
    }

    // run global stack registration to the template
    void ReconstructionModel::StackRegistrations3D(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations, RealImage& target){
        START_TIMING();

        InvertStackTransformations(stack_transformations);

        RealImage m_tmp = _evaluation_mask;
        TransformMask(target, m_tmp, RigidTransformation());
        target *= m_tmp;

        if (_verbose)
            SaveReconImageFile(target,"maskedTarget");
        

        RigidTransformation offset;
        ResetOrigin(target, offset);

        //register all stacks to the target
        Parallel::StackRegistrations p_reg(this, stacks, stack_transformations, target, offset);
        p_reg();

        for (int i=0; i<stack_transformations.size(); i++) {
            double tx, ty, tz, rx, ry, rz;
            tx = stack_transformations[i].GetTranslationX();
            ty = stack_transformations[i].GetTranslationY();
            tz = stack_transformations[i].GetTranslationZ();
            rx = stack_transformations[i].GetRotationX();
            ry = stack_transformations[i].GetRotationY();
            rz = stack_transformations[i].GetRotationZ();
            stringstream ss;
            ss << "- # " << i << " : " << tx << " " << ty << " " << tz << " | " << rx << " " << ry << " " << rz << endl;
            logger::debug(ss.str());
        }
        
        if(_debug)
            for (size_t i = 0; i < stacks.size(); i++)
                ReconstructionModel::SaveStackTransformation(stack_transformations[i], i);

        InvertStackTransformations(stack_transformations);


        END_TIMING("StackRegistrations");
    }

    //-------------------------------------------------------------------

    // restore the original slice intensities
    void ReconstructionModel::RestoreSliceIntensities() {
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            //calculate scaling factor
            const double& factor = _stack_factor[_stack_index[inputIndex]];//_average_value;

            // read the pointer to current slice
            RealPixel *p = _slices[inputIndex].Data();
            for (int i = 0; i < _slices[inputIndex].NumberOfVoxels(); i++)
                if (p[i] > 0)
                    p[i] /= factor;
        }
    }

    //-------------------------------------------------------------------

    double ReconstructionModel::CalculateSliceError(const int& index) {
        double sum = 0;
        if(((_slice_weight[index]>0.5) && _slice_inside[index])) {
            #pragma omp parallel for reduction(+: sum) collapse(3)
            for (int l=0; l<_slices[index].GetT(); l++) {
                for (int j=0; j<_slices[index].GetY(); j++) {
                    for (int i=0; i<_slices[index].GetX(); i++) {
                        if(_slices[index](i,j,0,l)>0 && _simulated_weights[index](i,j,0,l)>0.99)
                        {
                            double sliceValue = _slices[index](i,j,0,l) * exp(-(_bias[index])(i, j, 0)) * _scale[index];
                            double simValue = _simulated_slices[index](i,j,0,l);
                            double diff = sliceValue - simValue; 
                            sum += diff;
                        }     
                    }
                }
            }
        }
        return sum;
    }  


    double ReconstructionModel::Consistency(bool exclude_slices)
    {
        double sum = 0, num = 0, diff;
        
        #pragma omp parallel for reduction(+: diff, sum, num)
        for(int ind=0;ind<_slices.size();ind++) {
            if((!exclude_slices)||(_slice_weight[ind]>0.5)) {
                for (int i=0; i<_slices[ind].GetX();i++) {
                    for (int j=0; j<_slices[ind].GetY();j++) {
                        if(_slices[ind](i,j,0)>0) {
                            if(_simulated_weights[ind](i,j,0)>0.99){
                                diff = _slices[ind](i,j,0) - _simulated_slices[ind](i,j,0) * exp(_bias[ind](i, j, 0)) / _scale[ind];
                                sum+=diff*diff;
                                num++;
                            }
                        }
                    }
                }
            }
        }
        
        return sqrt(sum/num);
    }


    /***
     * Calculate Root Squared Error
     * as:
     * sqrt( sum ( (sliceVoxel - simSliceVoxel)^2 ) / num_voxels )
     */
    double ReconstructionModel::CalculateRSE(bool exclude_slices) {
        double sum = 0, num = 0;

        // slices
        for(int index=0; index < _slices.size(); index++) {
            if((!exclude_slices) || ((_slice_weight[index]>0.5) && _slice_inside[index])) {

                // voxels
                #pragma omp parallel for reduction(+: sum, num) collapse(3)
                for (int l=0; l<_slices[index].GetT(); l++) {
                    for (int j=0; j<_slices[index].GetY(); j++) {
                        for (int i=0; i<_slices[index].GetX(); i++) {
                            if(_slices[index](i,j,0,l)>0 && _simulated_weights[index](i,j,0,l)>0.99)
                            {
                                double sliceValue = _slices[index](i,j,0,l) * exp(-(_bias[index])(i, j, 0)) * _scale[index];
                                double simValue = _simulated_slices[index](i,j,0,l);
                                double diff = sliceValue - simValue; 
                                sum += diff*diff;
                                num++;
                            }
                        }
                    }
                }
            }
        }
        
        return sqrt(sum/num);
    }

    /**
     * Calculate Normalised Root Mean Squared Error
     * as:
     * sqrt( sum( (sliceVoxel - simSliceVoxel)^2 ) / num_voxels ) / 
     * ( sum(sliceVoxel) / num_voxels ) 
     */
    double ReconstructionModel::CalculateNRMSE(bool exclude_slices) {
        double sumDiff = 0, sumValue = 0, num = 0;

        // slices
        for(int index=0; index < _slices.size(); index++) {
            if((!exclude_slices) || ((_slice_weight[index]>0.5) && _slice_inside[index])) {

                // voxels
                #pragma omp parallel for reduction(+: sumDiff, sumValue, num) collapse(3)
                for (int l=0; l<_slices[index].GetT(); l++) {
                    for (int j=0; j<_slices[index].GetY(); j++) {
                        for (int i=0; i<_slices[index].GetX(); i++) {
                            if(_slices[index](i,j,0,l)>0 && _simulated_weights[index](i,j,0,l)>0.99)
                            {
                                double sliceValue = _slices[index](i,j,0,l) * exp(-(_bias[index])(i, j, 0)) * _scale[index];
                                double simValue = _simulated_slices[index](i,j,0,l);
                                double diff = sliceValue - simValue; 
                                sumValue += sliceValue;
                                sumDiff += diff*diff;
                                num++;
                            }
                        }
                    }
                }
            }
        }
        
        return sqrt(sumDiff/num) / (sumValue/num);
    }

    /**
     * Calculate Normalised Cross Coefficient
     * as:
     * sum((sliceVoxel - mean_sliceVoxel) * (simSliceVoxel - mean_simSliceVoxel)) /
     * sqrt( sum((sliceVoxel - mean_sliceVoxel)^2) * sum((simSliceVoxel - mean_simSliceVoxel)^2)) )
     */
    double ReconstructionModel::CalculateNCC(bool exclude_slices) {
        double sumSliceVal=0, countSliceVal=0, sumSimVal=0, countSimVal=0;
        double sliceAvg = 0, simAvg = 0;

        double diff_sum = 0;
        double sliceNormSq = 0;
        double simNormSq = 0;

        // Calculate means:

        // slices
        for(int index=0; index < _slices.size(); index++) {
            if((!exclude_slices) || ((_slice_weight[index]>0.5) && _slice_inside[index])) {

                // voxels
                #pragma omp parallel for reduction(+: sumSliceVal, countSliceVal, sumSimVal, countSimVal) collapse(3)
                for (int l=0; l<_slices[index].GetT(); l++) {
                    for (int j=0; j<_slices[index].GetY(); j++) {
                        for (int i=0; i<_slices[index].GetX(); i++) {
                            if(_slices[index](i,j,0,l)>0 && _simulated_weights[index](i,j,0,l)>0.99)
                            {
                                double sliceValue = _slices[index](i,j,0,l) * exp(-(_bias[index])(i, j, 0)) * _scale[index];
                                double simValue = _simulated_slices[index](i,j,0,l);

                                sumSliceVal += sliceValue;
                                countSliceVal++;
                                sumSimVal += simValue;
                                countSimVal++;
                            }
                        }
                    }
                }
            }
        }
        sliceAvg = sumSliceVal / countSliceVal;
        simAvg = sumSimVal / countSimVal;

        // Calculate NCC:

        // slices
        for(int index=0; index < _slices.size(); index++) {
            if((!exclude_slices) || ((_slice_weight[index]>0.5) && _slice_inside[index])) {

                // voxels
                #pragma omp parallel for reduction(+: diff_sum, sliceNormSq, simNormSq) collapse(3)
                for (int l=0; l<_slices[index].GetT(); l++) {
                    for (int j=0; j<_slices[index].GetY(); j++) {
                        for (int i=0; i<_slices[index].GetX(); i++) {
                            if(_slices[index](i,j,0,l)>0 && _simulated_weights[index](i,j,0,l)>0.99)
                            {
                                double sliceValue = _slices[index](i,j,0,l) * exp(-(_bias[index])(i, j, 0)) * _scale[index];
                                double simValue = _simulated_slices[index](i,j,0,l);

                                diff_sum += (sliceValue - sliceAvg) * (simValue - simAvg);
                                sliceNormSq += pow(sliceValue - sliceAvg, 2);
                                simNormSq += pow(simValue - simAvg, 2);

                            }
                        }
                    }
                }

            }
        }

        double ncc = 0;
        if (sliceNormSq * simNormSq > 0)
            ncc = diff_sum / sqrt(sliceNormSq * simNormSq);
        
        return ncc;
    }

    void ReconstructionModel::SaveMetaData(const boost::timer::cpu_timer& reconstruction_timer) {
        double ncc = CalculateNCC(true);
        double nrmse = CalculateNRMSE(true);
        double rse = CalculateRSE(true);
        int includedSlices = 0;
        int excludedSlices = 0;
        int outsideSlices = 0;
        SliceSummary(includedSlices, excludedSlices, outsideSlices);
        cout << "Creating json output" << endl;
        json j;
        j["evaluation"] = {
            {"ncc", ncc}, 
            {"nrmse", nrmse}, 
            {"rse", rse},
            {"count_slices_included", includedSlices},
            {"count_slices_excluded", excludedSlices},
            {"count_slices_outside_roi", outsideSlices},
            {"time_to_reconstruct", std::stof(reconstruction_timer.format(3,"%w"))}
        };
        std::map<int, std::vector<json>> stack_slice_info;
        for (int sliceIndex = 0; sliceIndex < _slices.size(); ++sliceIndex) {
            int stackIndex = _stack_index[sliceIndex];
            auto transform = _transformations[sliceIndex];
            json json_transform = json::object(
                                    {{"tx", transform.GetTranslationX() }, 
                                     {"ty", transform.GetTranslationY() }, 
                                     {"tz", transform.GetTranslationZ() }, 
                                     {"rx", transform.GetRotationX() }, 
                                     {"ry", transform.GetRotationY() }, 
                                     {"rz", transform.GetRotationZ() } });
            stack_slice_info[stackIndex].push_back( json::object({   
                {"id", sliceIndex},
                {"included", _slice_inside[sliceIndex] && _slice_weight[sliceIndex] >= 0.5 },
                {"excluded", _slice_inside[sliceIndex] && _slice_weight[sliceIndex] < 0.5 },
                {"outside_roi", !_slice_inside[sliceIndex] },
                {"stack", stackIndex},
                {"stack_position", _stack_loc_index[sliceIndex]},
                {"transformation", json_transform},
                {"error", CalculateSliceError(sliceIndex)},
                {"is_small", std::find(_small_slices.begin(), _small_slices.end(), sliceIndex) != _small_slices.end() }
            }));
        }
        
        std::vector<json> stacks;
        for (int stackIndex = 0; stackIndex < stack_slice_info.size(); ++stackIndex) {
            stacks.push_back(json::object({
                {"stack_number", stackIndex},
                {"slices", stack_slice_info[stackIndex]}
            }));
        }
        j["stacks"] = stacks;
        
        std::ofstream o("recon-metadata.json");
        o << std::setw(4) << j << std::endl;
    }


    //-------------------------------------------------------------------

    void ReconstructionModel::SliceSSIMMap(int inputIndex) {
        // Read slice
        RealImage slice = _slices[inputIndex];
        const RealImage& sim_slice = _simulated_slices[inputIndex];
        const RealImage& slice_mask = _slice_masks[inputIndex];

        // Adjust slice
        for (int l = 0; l < slice.GetT(); l++) {
            for (int j = 0; j < slice.GetY(); j++) {
                for (int i = 0; i < slice.GetX(); i++) {
                    if (InSliceMask(inputIndex, i, j)) {
                        slice(i, j, 0, l) *= exp(-(_bias[inputIndex](i, j, 0))) * _scale[inputIndex];
                    }
                }
            }
        }

        // Blur slice
        GaussianBlurring<RealPixel> gb(slice.GetXSize()*0.65);
        gb.Input(&slice);
        gb.Output(&slice);
        gb.Run();

        RealImage ssim_map = slice;
        ssim_map = 1;

        int shift = round(_local_SSIM_window_size/2);
        RealImage region_1, region_2, mask, region_mask;

        // For all slice voxels
        for (int t = 0; t < slice.GetT(); t++){
            for (int y = shift; y < (slice.GetY()-shift); y++) {
                for (int x = shift; x < (slice.GetX()-shift); x++) {
                    if (InSliceMask(inputIndex, x, y) && _structural_slice_weight[inputIndex] > 0) {
                        // Extract sub-region
                        region_1 = slice.GetRegion(x-shift, y-shift, 0, t, x+shift, y+shift, 1, t+1);
                        region_2 = sim_slice.GetRegion(x-shift, y-shift, 0, t, x+shift, y+shift, 1, t+1);
                        mask = slice_mask.GetRegion(x-shift, y-shift, 0, x+shift, y+shift, 1);

                        // Calculate region-mask
                        region_mask = region_1;
                        region_mask = 0;
                        for (int x=0; x<region_mask.GetX(); x++) {
                            for (int y=0; y<region_mask.GetY(); y++) {
                                int x2 = abs(shift-x)*abs(shift-x);
                                int y2 = abs(shift-y)*abs(shift-y);
                                int r = sqrt(x2+y2);
                                if (r < (shift+1))
                                    region_mask(x,y,0) = 1;
                            }
                        }
                        
                        // Combine mask and region mask
                        mask *= region_mask;

                        // Calculate local SSIM measure
                        double curSSIM = LocalSSIM(region_1, region_2, mask, t);

                        // Update SSIM map
                        ssim_map(x,y,0,t) = curSSIM;
                        if (ssim_map(x,y,0,t) < _local_SSIM_threshold)
                            ssim_map(x,y,0,t) = 0.1;
                        else
                            ssim_map(x,y,0,t) = 1;
                    } else {
                        ssim_map(x,y,0,t) = 1;
                    }
                }
            }
        }
        _slice_ssim_maps[inputIndex] = ssim_map;
        if (_verbose) SaveSliceFile(ssim_map, "ssim-map", inputIndex);

    }


    //-------------------------------------------------------------------

    // Implementation of SSIM between images
    double ReconstructionModel::LocalSSIM(const RealImage& slice, const RealImage& sim_slice, const RealImage& mask, const int t) {
        double C1 = 6.5025, C2 = 58.5225;
        double SSIM = 0, DSSIM = 0;
        
        double mu1=0, mu2=0, var1=0, var2=0, covar=0, num=0;
        double x_sq=0, y_sq=0, xy=0;                
        for (int yy=0; yy<slice.GetY(); yy++) {
            for (int xx=0; xx<slice.GetX(); xx++) {
                if (mask(xx,yy,0)==1 && slice(xx,yy,0)>0.01 && sim_slice(xx,yy,0)>0.01) {
                    mu1 += slice(xx, yy, 0, t);
                    mu2 += sim_slice(xx,yy,0,t);
                    num += 1;
                    x_sq += pow(slice(xx, yy, 0, t),2.0);
                    y_sq += pow(sim_slice(xx, yy, 0, t),2.0);
                    xy += slice(xx, yy, 0, t)*sim_slice(xx, yy, 0, t);
                 }
            }
        }
                
        mu1 = mu1/num;
        mu2 = mu2/num;
        var1 = (x_sq/num)-pow(mu1,2.0);
        var2 = (y_sq/num)-pow(mu2,2.0);
        covar = (xy/num)-mu1*mu2;
        double curSSIM = ((2*mu1*mu2+C1)*(2*covar+C2)) / ( (pow(mu1,2.)+pow(mu2,2.)+C1) * (var1+var2+C2) );
        
        return curSSIM;
    }

    
    //-------------------------------------------------------------------

    // run local structure-based oulier rejection step
    void ReconstructionModel::SStep() {
        START_TIMING();
        Parallel::SStep parallelSStep(this, _slices.size());
        parallelSStep();
        END_TIMING("SStep");
    }

    //-------------------------------------------------------------------

    // initialise / reset EM values
    void ReconstructionModel::InitializeEMValues() {
        START_TIMING();

        if (_no_masking_background) { //TODO - can this be removed?
            #pragma omp parallel for
            for (size_t i = 0; i < _slices.size(); i++) {
                //Initialise voxel weights and bias values
                const RealPixel *pi = _not_masked_slices[i].Data();
                RealPixel *pw = _weights[i].Data();
                RealPixel *pb = _bias[i].Data();
                for (int j = 0; j < _weights[i].NumberOfVoxels(); j++) {
                    pw[j] = pi[j] > -0.01 ? 1 : 0;
                    pb[j] = 0;
                }

                //Initialise slice weights
                _slice_weight[i] = 1;

                //Initialise scaling factors for intensity matching
                _scale[i] = 1;
            }
        } else {
            //#pragma omp parallel for
            for (size_t i = 0; i < _slices.size(); i++) {
                //Initialise voxel weights and bias values
                const RealPixel *pi = _slices[i].Data();
                RealPixel *pw = _weights[i].Data();
                RealPixel *pb = _bias[i].Data();
                for (int t = 0; t < _slices[i].GetT(); ++t)
                    for (int y = 0; y < _slices[i].GetY(); ++y)
                        for (int x = 0; x < _slices[i].GetX(); ++x) 
                            _weights[i](x,y,0,t) = InSliceMask(i,x,y) ? 1 : 0;
                
                for (int j = 0; j < _bias[i].NumberOfVoxels(); j++) {
                    pb[j] = 0;
                }

                //Initialise slice weights
                _slice_weight[i] = 1;

                //Initialise scaling factors for intensity matching
                _scale[i] = 1;
            }
        }

        //Force exclusion of slices predefined by user
        for (size_t i = 0; i < _force_excluded.size(); i++)
            if (_force_excluded[i] > 0 && _force_excluded[i] < _slices.size())
                _slice_weight[_force_excluded[i]] = 0;

        END_TIMING("InitializeEMValues");
    }

    //-------------------------------------------------------------------

    // run EStep for calculation of voxel-wise and slice-wise posteriors (weights)
    void ReconstructionModel::EStep() {
        Array<double> slice_potential(_slices.size());

        Parallel::EStep parallelEStep(this, slice_potential);
        parallelEStep();

        //To force-exclude slices predefined by a user, set their potentials to -1
        for (size_t i = 0; i < _force_excluded.size(); i++)
            if (_force_excluded[i] > 0 && _force_excluded[i] < _slices.size())
                slice_potential[_force_excluded[i]] = -1;

        //exclude slices identified as having small overlap with ROI, set their potentials to -1
        for (size_t i = 0; i < _small_slices.size(); i++)
            slice_potential[_small_slices[i]] = -1;

        //these are unrealistic scales pointing at misregistration - exclude the corresponding slices
        for (size_t inputIndex = 0; inputIndex < slice_potential.size(); inputIndex++)
            if (_scale[inputIndex] < 0.2 || _scale[inputIndex] > 5)
                slice_potential[inputIndex] = -1;

        //Calculation of slice-wise robust statistics parameters.
        //This is theoretically M-step,
        //but we want to use latest estimate of slice potentials
        //to update the parameters

        if (_verbose) {
            stringstream ss;
            ss << endl << "Slice potentials:";
            for (size_t inputIndex = 0; inputIndex < slice_potential.size(); inputIndex++)
                ss << " " << slice_potential[inputIndex];
            ss << endl;
            logger::trace(ss.str());
        }

        //Calculate means of the inlier and outlier potentials
        double sum = 0, den = 0, sum2 = 0, den2 = 0, maxs = 0, mins = 1;
        #pragma omp parallel for reduction(+: sum, sum2, den, den2) reduction(max: maxs) reduction(min: mins)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            if (slice_potential[inputIndex] >= 0) {
                //calculate means
                sum += slice_potential[inputIndex] * _slice_weight[inputIndex];
                den += _slice_weight[inputIndex];
                sum2 += slice_potential[inputIndex] * (1 - _slice_weight[inputIndex]);
                den2 += 1 - _slice_weight[inputIndex];

                //calculate min and max of potentials in case means need to be initalized
                if (slice_potential[inputIndex] > maxs)
                    maxs = slice_potential[inputIndex];
                if (slice_potential[inputIndex] < mins)
                    mins = slice_potential[inputIndex];
            }

        _mean_s = den > 0 ? sum / den : mins;
        _mean_s2 = den2 > 0 ? sum2 / den2 : (maxs + _mean_s) / 2;

        //Calculate the variances of the potentials
        sum = 0; den = 0; sum2 = 0; den2 = 0;
        #pragma omp parallel for reduction(+: sum, sum2, den, den2)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            if (slice_potential[inputIndex] >= 0) {
                sum += (slice_potential[inputIndex] - _mean_s) * (slice_potential[inputIndex] - _mean_s) * _slice_weight[inputIndex];
                den += _slice_weight[inputIndex];
                sum2 += (slice_potential[inputIndex] - _mean_s2) * (slice_potential[inputIndex] - _mean_s2) * (1 - _slice_weight[inputIndex]);
                den2 += 1 - _slice_weight[inputIndex];
            }

        //_sigma_s
        if (sum > 0 && den > 0) {
            //do not allow too small sigma
            _sigma_s = max(sum / den, _step * _step / 6.28);
        } else {
            _sigma_s = 0.025;
            if (_verbose) {
                if (sum <= 0)
                    logger::trace("All slices are equal. ");
                if (den < 0) //this should not happen
                    logger::trace("All slices are outliers. ");
                logger::trace((boost::format("Setting sigma to %1% ") % sqrt(_sigma_s)).str());
            }
        }

        //sigma_s2
        if (sum2 > 0 && den2 > 0) {
            //do not allow too small sigma
            _sigma_s2 = max(sum2 / den2, _step * _step / 6.28);
        } else {
            //do not allow too small sigma
            _sigma_s2 = max((_mean_s2 - _mean_s) * (_mean_s2 - _mean_s) / 4, _step * _step / 6.28);

            if (_verbose) {
                if (sum2 <= 0)
                    logger::trace("All slices are equal. ");
                if (den2 <= 0)
                    logger::trace("All slices are outliers. ");
                logger::trace((boost::format("Setting sigma_s2 to %1% ") % sqrt(_sigma_s2)).str());
            }
        }

        //Calculate slice weights
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            //Slice does not have any voxels in volumetric ROI
            if (slice_potential[inputIndex] == -1) {
                _slice_weight[inputIndex] = 0;
                continue;
            }

            //All slices are outliers or the means are not valid
            if (den <= 0 || _mean_s2 <= _mean_s) {
                _slice_weight[inputIndex] = 1;
                continue;
            }

            //likelihood for inliers
            const double gs1 = slice_potential[inputIndex] < _mean_s2 ? G(slice_potential[inputIndex] - _mean_s, _sigma_s) : 0;

            //likelihood for outliers
            const double gs2 = slice_potential[inputIndex] > _mean_s ? G(slice_potential[inputIndex] - _mean_s2, _sigma_s2) : 0;

            //calculate slice weight
            const double likelihood = gs1 * _mix_s + gs2 * (1 - _mix_s);
            if (likelihood > 0)
                _slice_weight[inputIndex] = gs1 * _mix_s / likelihood;
            else {
                if (slice_potential[inputIndex] <= _mean_s)
                    _slice_weight[inputIndex] = 1;
                if (slice_potential[inputIndex] >= _mean_s2)
                    _slice_weight[inputIndex] = 0;
                if (slice_potential[inputIndex] < _mean_s2 && slice_potential[inputIndex] > _mean_s) //should not happen
                    _slice_weight[inputIndex] = 1;
            }
        }

        //Update _mix_s this should also be part of MStep
        int num = 0; sum = 0;
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            if (slice_potential[inputIndex] >= 0) {
                sum += _slice_weight[inputIndex];
                num++;
            }

        if (num > 0)
            _mix_s = sum / num;
        else {
            cout << "All slices are outliers. Setting _mix_s to 0.9." << endl;
            _mix_s = 0.9;
        }

        if (_verbose) {
            stringstream ss;
            ss << setprecision(3);
            ss << "Slice robust statistics parameters: ";
            ss << "means: " << _mean_s << " " << _mean_s2 << "  ";
            ss << "sigmas: " << sqrt(_sigma_s) << " " << sqrt(_sigma_s2) << "  ";
            ss << "proportions: " << _mix_s << " " << 1 - _mix_s << endl;
            ss << "Slice weights:";
            for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
                ss << " " << _slice_weight[inputIndex];
            logger::trace(ss.str());
        }
    }

    //-------------------------------------------------------------------

    // run slice scaling
    void ReconstructionModel::Scale() {
        START_TIMING();

        Parallel::Scale parallelScale(this);
        parallelScale();

        if (_verbose) {
            stringstream ss;
            ss << setprecision(3);
            ss << "Slice scale =";
            for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
                ss << " " << _scale[inputIndex];
            logger::trace(ss.str());
        }

        END_TIMING("Scale");
    }

    //-------------------------------------------------------------------

    // run slice bias correction
    void ReconstructionModel::Bias() {
        START_TIMING();
        Parallel::Bias parallelBias(this);
        parallelBias();
        END_TIMING("Bias");
    }

    //-------------------------------------------------------------------

    // run MStep (RS)
    void ReconstructionModel::MStep(int iter) {
        Parallel::MStep parallelMStep(this);
        parallelMStep();
        const double sigma = parallelMStep.sigma;
        const double mix = parallelMStep.mix;
        const double num = parallelMStep.num;
        const double min = parallelMStep.min;
        const double max = parallelMStep.max;

        //Calculate sigma and mix
        if (mix > 0) {
            _sigma = sigma / mix;
        } else {
            throw runtime_error("Something went wrong: sigma=" + to_string(sigma) + " mix=" + to_string(mix));
        }
        if (_sigma < _step * _step / 6.28)
            _sigma = _step * _step / 6.28;
        if (iter > 1)
            _mix = mix / num;

        //Calculate m
        _m = 1 / (max - min);

        logger::trace((boost::format("Voxel-wise robust statistics parameters: sigma=%1% mix=%2% m=%3%") % sqrt(_sigma) % _mix % _m).str());
    }

    //-------------------------------------------------------------------

    // evaluation based on the number of excluded slices
    void ReconstructionModel::SliceSummary(int& included_count, int& excluded_count, int& outside_count) {
        included_count = 0, excluded_count = 0, outside_count = 0;

        string included, excluded, outside;
        for (size_t i = 0; i < _slices.size(); i++) {
            if (_slice_inside[i]) {
                if (_slice_weight[i] >= 0.5) {
                    included += " " + to_string(i);
                    included_count++;
                } else {
                    excluded += " " + to_string(i);
                    excluded_count++;
                }
            } else {
                outside += " " + to_string(i);
                outside_count++;
            }
        }
        if (_debug) {
            logger::info("Included slices: " + included);
            logger::info("Excluded slices: " + excluded);
            logger::info("Outside slices: " + outside);
        }
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::ReadTransformations(const char *folder, size_t file_count, Array<RigidTransformation>& transformations) {
        if (_slices.empty())
            throw runtime_error("Please create slices before reading transformations!");

        ClearAndReserve(transformations, file_count);
        for (size_t i = 0; i < file_count; i++) {
            const string path = (boost::format("%1%/transformation%2%.dof") % (folder ? folder : ".") % i).str();
            Transformation *transformation = Transformation::New(path.c_str());
            unique_ptr<RigidTransformation> rigidTransf(dynamic_cast<RigidTransformation*>(transformation));
            transformations.push_back(*rigidTransf);
            cout << path << endl;
        }
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::ReadTransformations(const char *folder) {
        cout << "Reading transformations:" << endl;
        ReadTransformations(folder, _slices.size(), _transformations);
    }

    //-------------------------------------------------------------------
    

    /// ***********************
    /// Save functions
    /// ***********************


    /// Output as stack Nifti file
    void ReconstructionModel::SaveStackFile(const RealImage& image, 
                        const string& suffix, 
                        const int stackNum ,
                        const int iterNum ,
                        const int recIterNum) {
        logger::trace("Saving " + suffix + " iter=" + to_string(iterNum) + " recIterNum=" + to_string(recIterNum));
        string filename;
        if (stackNum >= 0 && iterNum >= 0 && recIterNum >= 0)
            filename = (boost::format("stack%1%-it%2%-sr%3%-%4%.nii.gz") % stackNum % iterNum % recIterNum % suffix).str();
        else if  (stackNum >= 0 &&  iterNum >= 0 )
            filename = (boost::format("stack%1%-it%2%-%3%.nii.gz") % stackNum % iterNum % suffix).str();
        else if (stackNum >= 0)
            filename = (boost::format("stack%1%-%2%.nii.gz") % stackNum % suffix).str();
        else
            filename = (boost::format("stack-%1%.nii.gz") % suffix).str();
        logger::trace("Saving " + filename);
        image.Write(filename.c_str());
        
    }

    /// Output image as nifti file 
    void ReconstructionModel::SaveReconImageFile(const RealImage& image,
                                        const string& suffix, 
                                        const int iterNum,
                                        const int recIterNum) {
        string filename;
        if (iterNum >= 0 && recIterNum >= 0)
            filename = (boost::format("reconstruction-it%1%-sr%2%-%3%.nii.gz") % iterNum % recIterNum % suffix).str();
        else if  (iterNum >= 0 )
            filename = (boost::format("reconstruction-it%1%-%2%.nii.gz") % iterNum % suffix).str();
        else
            filename = (boost::format("reconstruction-%1%.nii.gz") % suffix).str();
        logger::trace("Saving " + filename);
        image.Write(filename.c_str());
    }
    
    void ReconstructionModel::SaveReconImageFile(const RealImage& image,
                                    const int iterNum,
                                    const int recIterNum) {
        string filename;
        if (recIterNum >= 0)
            filename = (boost::format("reconstruction-it%1%-sr%2%.nii.gz") % iterNum % recIterNum).str();
        else
            filename = (boost::format("reconstruction-it%1%.nii.gz") % iterNum).str();
        logger::trace("Saving " + filename);
        image.Write(filename.c_str());
    }


    void ReconstructionModel::SaveSliceFile(const RealImage& image, const string& suffix, const int sliceNum) {
        string filename;
        filename = (boost::format("slice%1%-%2%.nii.gz") % sliceNum % suffix).str();
        logger::trace("Saving " + filename);
        image.Write(filename.c_str());
    }

    void ReconstructionModel::SaveStacks(const Array<RealImage>& save_stacks, const string& name) {
        logger::trace("Saving stacks as "+name);
        #pragma omp parallel for
        for (size_t i = 0; i < save_stacks.size(); i++)
            SaveStackFile(save_stacks[i],name,i);
    }

    /// Output stack transform as dof file
    void ReconstructionModel::SaveStackTransformation(const RigidTransformation& transform, const int stackNum) {
        transform.Write((boost::format("stack%1%-transformation.dof") % stackNum ).str().c_str());
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::SaveSlices() {
        logger::trace("Saving slices");
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            SaveSliceFile(_slices[inputIndex], "slice", inputIndex);
    }
    
    //-------------------------------------------------------------------

    void ReconstructionModel::SaveSimulatedSlices() {
        logger::trace("Saving simulated slices ... ");
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            SaveSliceFile(_simulated_slices[inputIndex], "simulated", inputIndex);
    }


    //-------------------------------------------------------------------

    void ReconstructionModel::SaveSliceMasks() {
        logger::trace("Saving slice masks slices ... ");
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            SaveSliceFile(_slice_masks[inputIndex], "mask", inputIndex);
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::SaveWeights() {
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            SaveSliceFile(_weights[inputIndex], "weights", inputIndex);
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::SaveTransformations() {
        #pragma omp parallel for
        for (size_t i = 0; i < _transformations.size(); i++)
            _transformations[i].Write((boost::format("slice%1%-transformation.dof") % i).str().c_str());
    }

    //-------------------------------------------------------------------

    void ReconstructionModel::SaveProbabilityMap(int i) {
        _brain_probability.Write((boost::format("slice%1%-probabilitymap.nii") % i).str().c_str());
    }

    void ReconstructionModel::GetSimulatedSlicesAsStacks(const vector<RealImage>& stacks, vector<RealImage>& simulatedStacks){
        logger::trace("GetSimulatedSlicesAsStacks ");

        simulatedStacks.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++)
            simulatedStacks.push_back(RealImage(stacks[i].Attributes()));

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _simulated_slices.size(); inputIndex++)
            for (int t = 0; t < _simulated_slices[inputIndex].GetT(); t++)
                for (int j = 0; j < _simulated_slices[inputIndex].GetY(); j++)
                    for (int i = 0; i < _simulated_slices[inputIndex].GetX(); i++)
                        simulatedStacks[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], t) = _simulated_slices[inputIndex](i, j, 0, t);
    }

    
    void ReconstructionModel::GetScaledSimulatedSlicesAsStacks(const vector<RealImage>& stacks, vector<RealImage>& simulatedStacks){
        logger::trace("GetSimulatedSlicesAsStacks ");

        simulatedStacks.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++)
            simulatedStacks.push_back(RealImage(stacks[i].Attributes()));

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _simulated_slices.size(); inputIndex++)
            for (int t = 0; t < _simulated_slices[inputIndex].GetT(); t++)
                for (int j = 0; j < _simulated_slices[inputIndex].GetY(); j++)
                    for (int i = 0; i < _simulated_slices[inputIndex].GetX(); i++)
                        simulatedStacks[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], t) = 
                            _simulated_slices[inputIndex](i, j, 0, t)  * exp(_bias[inputIndex](i, j, 0)) / _scale[inputIndex];
    }

    //-------------------------------------------------------------------


    /// Save simulated slices
    void ReconstructionModel::SaveSimulatedSlicesAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        logger::trace("Saving simulated slices as stack for iter " + to_string(iter));
        SaveAsStack(_simulated_slices, stacks, "simulated", iter, rec_iter);
    }
    
    /// Save simulated slices
    void ReconstructionModel::SaveScaledSimulatedSlicesAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        vector<RealImage> scaledSimulatedStacks;
        GetScaledSimulatedSlicesAsStacks(stacks, scaledSimulatedStacks);
        logger::trace("Saving simulated slices as stack for iter " + to_string(iter));
        //SaveAsStack(scaledSimulatedStacks, stacks, "scaledSimulated", iter, rec_iter);
        
        for (size_t i = 0; i < stacks.size(); i++) {
            string filename;
                SaveStackFile(scaledSimulatedStacks[i],"scaledSimulated",i,iter,rec_iter);
        }
    }

    /// Save simulated weights
    void ReconstructionModel::SaveSimulatedWeights(const Array<RealImage>& stacks, int iter, int rec_iter) {
        SaveAsStack(_simulated_weights, stacks, "simweightstack", iter, rec_iter);
    }
    
    void ReconstructionModel::SaveSlicesAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        SaveAsStack(_slices, stacks, "slices", iter, rec_iter);
    }

    // Save weights
    void ReconstructionModel::SaveWeightsAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        SaveAsStack(_weights, stacks, "weight", iter, rec_iter);
    }


    void ReconstructionModel::SaveSSIMMapsAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
         SaveAsStack(_slice_ssim_maps, stacks, "ssim_maps", iter, rec_iter);
    }
    
    // mask slices based on the reconstruction mask
    void ReconstructionModel::MaskSlices() {
        //Check whether we have a mask
        if (!_have_mask) {
            cerr << "Could not mask slices because no mask has been set." << endl;
            return;
        }

        if (_no_masking_background)
            _not_masked_slices.insert(_not_masked_slices.end(), _slices.begin(), _slices.end());

        Utility::MaskSlices(_slices, _reconstructed_mask, [&](size_t index, double& x, double& y, double& z) {
            _transformations[index].Transform(x, y, z);
        });
    }


    void ReconstructionModel::CreateSliceMasks(){
        START_TIMING();
        
        //Check whether we have a mask
        if (!_have_mask) {
            cerr << "Could not mask slices because no mask has been set." << endl;
            return;
        }

        Utility::CreateSliceMasks(_slices, _slice_masks, _reconstructed_mask, 
            [&](size_t index, double& x, double& y, double& z) { _transformations[index].Transform(x, y, z); }
            );
            
        END_TIMING("CreateSliceMasks");
    }

    //-------------------------------------------------------------------


    // -----------------------------------------------------------------------------
    // Determine Reconstructed Spatial Resolution
    // -----------------------------------------------------------------------------
    //determine resolution of volume to reconstruct
    double ReconstructionModel::GetReconstructedResolutionFromTemplateStack(const RealImage& stack) {
        double dx, dy, dz, d;
        stack.GetPixelSize(&dx, &dy, &dz);
        if (dx <= dy && dx <= dz)
            d = dx;
        else if (dy <= dz)
            d = dy;
        else
            d = dz;
        return d;
    }

    // -----------------------------------------------------------------------------
    // Get Slice-Location Transformations
    // -----------------------------------------------------------------------------
    void ReconstructionModel::ReadSliceTransformations(const char *folder) {
        stringstream ss;
        ss<< "Reading slice transformations from: " << folder;
        logger::info(ss.str());
        Array<RigidTransformation> loc_transformations;
        ReadTransformations(folder, _loc_index.back() + 1, loc_transformations);

        // Assign transformations to single-frame images
        ClearAndReserve(_transformations, _slices.size());
        for (size_t i = 0; i < _slices.size(); i++)
            _transformations.push_back(loc_transformations[_loc_index[i]]);
    }

    // -----------------------------------------------------------------------------
    // Initialise Reconstructed Volume from Static Mask
    // -----------------------------------------------------------------------------
    // Create zero image as a template for reconstructed volume
    void ReconstructionModel::CreateTemplateFromStaticMask(const RealImage& mask, 
                                                                      double resolution,
                                                                      int temporalSize,
                                                                      double temporalResolution) {
        // Get mask attributes - image size and voxel size
        ImageAttributes attr = mask.Attributes();

        // Set temporal dimension
        attr._t = temporalSize;
        attr._dt = temporalResolution;

        // Create volume and initialise volume values
        RealImage volume4D(attr);

        // Resample to specified resolution
        unique_ptr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(Interpolation_NN));

        Resampling<RealPixel> resampling(resolution, resolution, resolution);
        resampling.Input(&volume4D);
        resampling.Output(&volume4D);
        resampling.Interpolator(interpolator.get());
        resampling.Run();

        // Set reconstructed 4D volume
        _reconstructed = move(volume4D);
        _template_created = true;

        logger::debug("CreateTemplateFromStaticMask created template 4D volume with "
            + to_string(volume4D.GetT()) + " time points and " + to_string(volume4D.TSize()) +  "s temporal resolution.");
    }

    void ReconstructionModel::CreateTemplate(ImageAttributes attr, 
                                             const RealImage& templateImage,
                                             double resolution,
                                             int temporalSize,
                                             double temporalResolution) {

        // Set temporal dimension
        attr._t = temporalSize;
        attr._dt = temporalResolution;

        // Create volume and initialise volume values
        RealImage volume4D(attr);

        // Calculate padding requirements
        RealPixel smin, smax;
        templateImage.GetMinMax(&smin, &smax);
        int pad;
        bool usePad = false;
        if (smin < -0.1) {
            pad = -1;
            usePad = true;
        } else if (smin < 0.1) {
            pad = 0;
            usePad = true;
        }

        // Resample
        if (usePad) {
            GenericLinearInterpolateImageFunction<RealImage> interpolator;
            ResamplingWithPadding<RealPixel> resampler(resolution, resolution, resolution, pad);
            resampler.Input(&templateImage);
            resampler.Output(&volume4D);
            resampler.Interpolator(&interpolator);
            resampler.Run();
        }
        else {
            unique_ptr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(Interpolation_Linear));
            Resampling<RealPixel> resampler(resolution, resolution, resolution);
            resampler.Input(&templateImage);
            resampler.Output(&volume4D);
            resampler.Interpolator(interpolator.get());
            resampler.Run();
        }

        // Set reconstructed 4D volume
        _reconstructed = move(volume4D);
        _template_created = true;
        
        if (_verbose) SaveReconImageFile(_reconstructed,"template");
    }

    // -----------------------------------------------------------------------------
    // Match Stack Intensities With Masking
    // -----------------------------------------------------------------------------
    void ReconstructionModel::MatchStackIntensitiesWithMasking(Array<RealImage>& stacks, const Array<RigidTransformation>& stack_transformations, double averageValue, bool together) {
        logger::debug("Matching intensities of stacks. ");

        cout << setprecision(6);

        //Calculate the averages of intensities for all stacks
        RealImage m;
        Array<double> stack_average;
        stack_average.reserve(stacks.size());

        //remember the set average value
        _average_value = averageValue;

        //averages need to be calculated only in ROI
        for (size_t ind = 0; ind < stacks.size(); ind++) {
            ImageAttributes attr = stacks[ind].Attributes();
            attr._t = 1;
            m.Initialize(attr);
            double sum = 0, num = 0;

            #pragma omp parallel for reduction(+: sum, num)
            for (int i = 0; i < stacks[ind].GetX(); i++)
                for (int j = 0; j < stacks[ind].GetY(); j++)
                    for (int k = 0; k < stacks[ind].GetZ(); k++) {
                        //image coordinates of the stack voxel
                        double x = i;
                        double y = j;
                        double z = k;
                        //change to world coordinates
                        stacks[ind].ImageToWorld(x, y, z);
                        //transform to template (and also _mask) space
                        stack_transformations[ind].Transform(x, y, z);
                        //change to mask image coordinates - mask is aligned with template
                        _reconstructed_mask.WorldToImage(x, y, z);
                        x = round(x);
                        y = round(y);
                        z = round(z);
                        //if the voxel is inside mask ROI include it
                        if (_reconstructed_mask.IsInside(x, y, z)) {
                            if (_reconstructed_mask(x, y, z) == 1) {
                                m(i, j, k) = 1;
                                for (int f = 0; f < stacks[ind].GetT(); f++) {
                                    sum += stacks[ind](i, j, k, f);
                                    num++;
                                }
                            }
                        }
                    }
            if (_verbose)
                SaveStackFile(m,"maskForIntensityMatching",ind); 

            //calculate average for the stack
            if (num > 0) {
                stack_average.push_back(sum / num);
            } else {
                throw runtime_error("Stack " + to_string(ind) + " has no overlap with ROI");
            }
        }

        double global_average = 0;
        if (together) {
            for (size_t i = 0; i < stack_average.size(); i++)
                global_average += stack_average[i];
            global_average /= stack_average.size();
        }

        if (_debug) {
            std::stringstream ss;
            ss << setprecision(6);
            ss << "Stack average intensities are ";
            for (size_t ind = 0; ind < stack_average.size(); ind++)
                ss << stack_average[ind] << " ";
            ss << "\nThe new average value is " << averageValue;
            logger::debug(ss.str());
        }

        //Rescale stacks
        ClearAndReserve(_stack_factor, stacks.size());
        for (size_t ind = 0; ind < stacks.size(); ind++) {
            double factor = averageValue / (together ? global_average : stack_average[ind]);
            _stack_factor.push_back(factor);

            RealPixel *ptr = stacks[ind].Data();
            #pragma omp parallel for
            for (int i = 0; i < stacks[ind].NumberOfVoxels(); i++) {
                if (ptr[i] > 0)
                    ptr[i] *= factor;
            }
        }

        if (_verbose) {
            #pragma omp parallel for
            for (size_t ind = 0; ind < stacks.size(); ind++)
                SaveStackFile(stacks[ind],"rescaled", ind);

            stringstream ss;
            ss << "Slice intensity factors are ";
            for (size_t ind = 0; ind < _stack_factor.size(); ind++)
                ss << _stack_factor[ind] << " ";
            ss << "\nThe new average value is " << averageValue;
            logger::debug(ss.str());
        }

        cout << setprecision(3);
    }

    // -----------------------------------------------------------------------------
    // Create Slices and Associated Transformations
    // -----------------------------------------------------------------------------
    void ReconstructionModel::CreateSlicesAndTransformations(const Array<RealImage>& stacks, const Array<RigidTransformation>& stack_transformations, const Array<double>& thickness, const Array<RealImage>& probability_maps) {
        logger::debug("CreateSlicesAndTransformations");
        START_TIMING();

        // Reset and allocate memory
        size_t reserve_size, reserve_size_3d;
        reserve_size = stacks.size() * stacks[0].Attributes()._z * stacks[0].Attributes()._t;
        reserve_size_3d = stacks.size() * stacks[0].Attributes()._z;


        ClearAndReserve(_slices, reserve_size);
        ClearAndReserve(_simulated_slices, reserve_size);
        ClearAndReserve(_simulated_weights, reserve_size);
        ClearAndReserve(_simulated_inside, reserve_size);
        ClearAndReserve(_zero_slices, reserve_size);
        ClearAndReserve(_stack_index, reserve_size);
        ClearAndReserve(_loc_index, reserve_size);
        ClearAndReserve(_stack_loc_index, reserve_size);
        ClearAndReserve(_transformations, reserve_size);
        Clear(_slice_excluded, reserve_size);
        ClearAndReserve(_slice_masks, reserve_size);
        ClearAndReserve(_structural_slice_weight, reserve_size);
        ClearAndReserve(_slice_ssim_maps, reserve_size);

        if (!probability_maps.empty())
            ClearAndReserve(_probability_maps, reserve_size);

        for (size_t i = 0, loc_index = 0; i < stacks.size(); i++) {
            //image attributes contain image and voxel size
            const ImageAttributes& attr = stacks[i].Attributes();

            //attr._z is number of slices in the stack
            for (int j = 0; j < attr._z; j++, loc_index++) {
                //create slice by selecting the appropriate region of the stack
                RealImage slice = stacks[i].GetRegion(0, 0, j, 0, attr._x, attr._y, j + 1, attr._t);
                //set correct voxel size in the stack. Z size is equal to slice thickness.
                slice.PutPixelSize(attr._dx, attr._dy, thickness[i], attr._dt);
                //remember the slice
                _slices.push_back(slice);
                _slice_masks.push_back(slice.GetRegion(0,0,0,0,attr._x,attr._y,1,1) = 1); // Default mask with all voxels set to 1.
                _simulated_slices.push_back(slice);
                _simulated_weights.push_back(slice);
                _simulated_inside.push_back(move(slice));
                _slice_ssim_maps.push_back(slice);
                _zero_slices.push_back(1);
                //remember stack indices for this slice
                _stack_index.push_back(i);
                _loc_index.push_back(loc_index);
                _stack_loc_index.push_back(j);
                //initialize slice transformation with the stack transformation
                _transformations.push_back(stack_transformations[i]);
                //initialise slice exclusion flags
                _slice_excluded.push_back(0);
                if (!probability_maps.empty()) {
                    RealImage proba = probability_maps[i].GetRegion(0, 0, j, 0, attr._x, attr._y, j + 1, attr._t);
                    proba.PutPixelSize(attr._dx, attr._dy, thickness[i], attr._dt);
                    _probability_maps.push_back(move(proba));
                }   
                _structural_slice_weight.push_back(1);
            }
        }

        logger::debug("Number of images: " + to_string(_slices.size()));

        //set excluded slices
        for (size_t i = 0; i < _force_excluded.size(); i++)
            _slice_excluded[_force_excluded[i]] = true;

        //set excluded stacks
        for (size_t i = 0; i < _force_excluded_stacks.size(); i++)
            for (size_t j = 0; j < _slices.size(); j++)
                if (_force_excluded_stacks[i] == _stack_index[j])
                    _slice_excluded[j] = true;
        
        END_TIMING("CreateSlicesAndTransformations");
    }

    // -----------------------------------------------------------------------------
    // Calculate Transformation Matrix Between Slices and Voxels
    // -----------------------------------------------------------------------------
    void ReconstructionModel::CoeffInit() {
        START_TIMING();

        //resize slice-volume matrix from previous iteration
        ClearAndResize(_volcoeffs, _slices.size());

        //resize indicator of slice having and overlap with volumetric mask
        ClearAndResize(_slice_inside, _slices.size());

        logger::trace("Initialising matrix coefficients... ");
        Parallel::CoeffInit coeffinit(this);
        coeffinit();
        logger::trace( "done.");

        //prepare image for volume weights, will be needed for Gaussian Reconstruction
        logger::trace("Computing 4D volume weights...");
        _volume_weights.Initialize(_reconstructed.Attributes());

        logger::trace("    ... for input slice: ");
        // Do not parallelise: It would cause data inconsistencies
        stringstream ss;
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            if (_verbose)
                ss << inputIndex << ", ";

            bool excluded = false;

            for (size_t fe = 0; fe < _force_excluded.size(); fe++) {
                if (inputIndex == _force_excluded[fe]) {
                    excluded = true;
                    break;
                }
            }

            if (_structural_slice_weight[inputIndex] < 0.5)
                excluded = true;

            if (!excluded) {
            // Do not parallelise: It would cause data inconsistencies
            for (size_t i = 0; i < _volcoeffs[inputIndex].size(); i++)
                for (size_t j = 0; j < _volcoeffs[inputIndex][i].size(); j++)
                    for (size_t k = 0; k < _volcoeffs[inputIndex][i][j].size(); k++) {
                        const POINT3D& p = _volcoeffs[inputIndex][i][j][k];
                        // Future potential improvement: should _volume_weight be 3D only if repeating the data?
                        for (int t = 0; t < _reconstructed.GetT(); t++)
                            _volume_weights(p.x, p.y, p.z, t) += p.value;
                    }
            }
        }
        if (_verbose)
            ss << "\b\b" << endl;
        logger::trace(ss.str());

        //find average volume weight to modify alpha parameters accordingly
        double sum = 0;
        int num = 0;
        #pragma omp parallel for reduction(+: sum, num)
        for (int i = 0; i < _volume_weights.GetX(); i++)
            for (int j = 0; j < _volume_weights.GetY(); j++)
                for (int k = 0; k < _volume_weights.GetZ(); k++)
                    if (_reconstructed_mask(i, j, k) == 1)
                        for (int f = 0; f < _volume_weights.GetT(); f++) {
                            sum += _volume_weights(i, j, k, f);
                            num++;
                        }

        _average_volume_weight = sum / num;

        logger::trace((boost::format("Average volume weight is  %1%") % _average_volume_weight).str());

        END_TIMING("CoeffInit");
    }

    // -----------------------------------------------------------------------------
    // PSF-Weighted Reconstruction
    // -----------------------------------------------------------------------------
    void ReconstructionModel::GaussianReconstruction() {
        START_TIMING();
        // SaveReconImageFile(_volume_weights,"gr-volumeWeights");

        logger::trace("Gaussian reconstruction ... ");

        RealImage slice;
        Array<int> voxel_num;
        voxel_num.reserve(_slices.size() * _slices[0].GetX() * _slices[0].GetY() * _slices[0].GetT());


        //clear _reconstructed image
        memset(_reconstructed.Data(), 0, sizeof(RealPixel) * _reconstructed.NumberOfVoxels());

        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            if (_slice_excluded[inputIndex])
                continue;

            //copy the current slice
            slice = _slices[inputIndex];
            const RealImage& mask = _slice_masks[inputIndex];

            //alias the current bias image
            const RealImage& b = _bias[inputIndex];
            //read current scale factor
            const double scale = _scale[inputIndex];

            int slice_vox_num = 0;

            //Distribute slice intensities to the volume
            for (int t = 0; t < slice.GetT(); t++) {
                for (size_t i = 0; i < _volcoeffs[inputIndex].size(); i++) {
                    for (size_t j = 0; j < _volcoeffs[inputIndex][i].size(); j++) {
                        if (InSliceMask(inputIndex,i, j)) {
                            //biascorrect and scale the slice
                            slice(i, j, 0, t) *= exp(-b(i, j, 0)) * scale;

                            //number of volume voxels with non-zero coefficients
                            //for current slice voxel
                            const size_t n = _volcoeffs[inputIndex][i][j].size();

                            //calculate num of vox in a slice that have overlap with roi
                            if (n > 0)
                                slice_vox_num++;

                            //add contribution of current slice voxel to all voxel volumes
                            //to which it contributes
                            for (size_t k = 0; k < n; k++) {
                                const POINT3D& p = _volcoeffs[inputIndex][i][j][k];
                                _reconstructed(p.x, p.y, p.z, t) += p.value * slice(i, j, 0, t);
                            }
                        }
                    }
                }
            }

            voxel_num.push_back(slice_vox_num);
        } //end of loop for a slice inputIndex

        //normalize the volume by proportion of contributing slice voxels
        //for each volume voxel
        _reconstructed /= _volume_weights;

        logger::trace("Gaussian reconstruction done.");

        // find slices with small overlap with ROI and exclude them.
        //find median
        Array<int> voxel_num_tmp = voxel_num;
        int median = round(voxel_num_tmp.size() * 0.5) - 1;
        nth_element(voxel_num_tmp.begin(), voxel_num_tmp.begin() + median, voxel_num_tmp.end());
        median = voxel_num_tmp[median];

        //remember slices with small overlap with ROI
        ClearAndReserve(_small_slices, voxel_num.size());
        for (size_t i = 0; i < voxel_num.size(); i++)
            if (voxel_num[i] < 0.1 * median)
                _small_slices.push_back(i);

        if (_verbose) {
            stringstream ss;
            ss << "Small slices:";
            for (size_t i = 0; i < _small_slices.size(); i++)
                ss << " " << _small_slices[i];
            ss << endl;
            logger::trace(ss.str());
        }

        END_TIMING("GaussianReconstruction");
    }

    // -----------------------------------------------------------------------------
    // Calculate Corrected Slices
    // -----------------------------------------------------------------------------
    void ReconstructionModel::CalculateCorrectedSlices() {
        logger::trace("Calculating corrected slices ...");

        Parallel::CalculateCorrectedSlices parallelCalculateCorrectedSlices(this);
        parallelCalculateCorrectedSlices();

        logger::trace("done.");
    }

    // -----------------------------------------------------------------------------
    // Simulate Slices
    // -----------------------------------------------------------------------------
    void ReconstructionModel::SimulateSlices() {
        START_TIMING();

        logger::trace("Simulating slices ... ");

        Parallel::SimulateSlices parallelSimulateSlices(this);
        parallelSimulateSlices();

        logger::trace("done.");

        END_TIMING("SimulateSlices");
    }

    // -----------------------------------------------------------------------------
    // Simulate Stacks
    // -----------------------------------------------------------------------------
    void ReconstructionModel::SimulateStacks() {
        logger::trace("Simulating stacks ... ");

        //Initialise simulated images
        #pragma omp parallel for
        for (size_t i = 0; i < _slices.size(); i++)
            memset(_simulated_slices[i].Data(), 0, sizeof(RealPixel) * _simulated_slices[i].NumberOfVoxels());

        //Initialise indicator of images having overlap with volume
        _slice_inside = deque<bool>(_slices.size(), true);

        //Simulate images
        Parallel::SimulateStacks simulatestacks(this);
        simulatestacks();

        logger::trace("done.");
    }

    // -----------------------------------------------------------------------------
    // Calculate Error
    // -----------------------------------------------------------------------------
    void ReconstructionModel::CalculateError() {
        logger::trace("Calculating error ... ");

        Parallel::CalculateError parallelCalculateError(this);
        parallelCalculateError();

        logger::trace("done.");
    }

    // -----------------------------------------------------------------------------
    // Normalise Bias
    // -----------------------------------------------------------------------------
    void ReconstructionModel::NormaliseBias(int iter, int rec_iter) {
        START_TIMING();

        logger::trace("Normalise Bias ... ");

        Parallel::NormaliseBias parallelNormaliseBias(this);
        parallelNormaliseBias();
        RealImage& bias = parallelNormaliseBias.bias;
        RealImage& volweight3d = parallelNormaliseBias.volweight3d;

        // normalize the volume by proportion of contributing slice voxels for each volume voxel
        bias /= volweight3d;

        logger::trace("done.");

        StaticMaskVolume(bias, _reconstructed_mask, 0);
        RealImage m = _reconstructed_mask;
        GaussianBlurring<RealPixel> gb(_sigma_bias);
        gb.Input(&bias);
        gb.Output(&bias);
        gb.Run();
        gb.Input(&m);
        gb.Output(&m);
        gb.Run();

        bias /= m;

        if (_verbose)
            SaveReconImageFile(bias,"averageBias",iter,rec_iter);

        #pragma omp parallel for
        for (int i = 0; i < _reconstructed.GetX(); i++)
            for (int j = 0; j < _reconstructed.GetY(); j++)
                for (int k = 0; k < _reconstructed.GetZ(); k++)
                    for (int f = 0; f < _reconstructed.GetT(); f++)
                        if (_reconstructed(i, j, k, f) != -1)
                            _reconstructed(i, j, k, f) /= exp(-bias(i, j, k));

        END_TIMING("NormaliseBias");
    }

    // -----------------------------------------------------------------------------
    // Slice-to-Volume Registration
    // -----------------------------------------------------------------------------
    void ReconstructionModel::SliceToVolumeRegistration(int iter, bool maskSlices) {
        START_TIMING();
        logger::trace("SliceToVolumeRegistration");

        Parallel::SliceToVolumeRegistration registration(this, iter, maskSlices);
        registration();

        END_TIMING("SliceToVolumeRegistration");
    }

    // -----------------------------------------------------------------------------
    // Volume-to-Volume Registration
    // -----------------------------------------------------------------------------
    void ReconstructionModel::VolumeToVolumeRegistration(const GreyImage& target, const GreyImage& source, RigidTransformation& rigidTransf) {
        ParameterList params;
        Insert(params, "Transformation model", "Rigid");
        Insert(params, "Image (dis-)similarity measure", "NMI");
        if (_nmi_bins > 0)
            Insert(params, "No. of bins", _nmi_bins);
        Insert(params, "Background value for image 1", -1);
        Insert(params, "Background value for image 2", -1);

        GenericRegistrationFilter registration;
        registration.Parameter(params);
        registration.Input(&target, &source);
        Transformation *dofout = nullptr;
        registration.Output(&dofout);
        registration.InitialGuess(&rigidTransf);
        registration.GuessParameter();
        registration.Run();

        unique_ptr<RigidTransformation> rigidTransf_dofout(dynamic_cast<RigidTransformation*>(dofout));
        rigidTransf = *rigidTransf_dofout;
    }

    // -----------------------------------------------------------------------------
    // Common Calculate Function
    // -----------------------------------------------------------------------------
    double ReconstructionModel::Calculate(const RigidTransformation& drift, double (ReconstructionModel::*Calculate)()) {
        //Initialise
        Matrix d = drift.GetMatrix();

        //Remove drift
        #pragma omp parallel for
        for (size_t i = 0; i < _transformations.size(); i++)
            _transformations[i].PutMatrix(d * _transformations[i].GetMatrix());

        //Calculate
        const double mean = (this->*Calculate)();

        //Return drift
        d.Invert();
        #pragma omp parallel for
        for (size_t i = 0; i < _transformations.size(); i++)
            _transformations[i].PutMatrix(d * _transformations[i].GetMatrix());

        //Output
        return mean;
    }

    // -----------------------------------------------------------------------------
    // Calculate Displacement
    // -----------------------------------------------------------------------------
    double ReconstructionModel::CalculateDisplacement() {
        logger::debug("CalculateDisplacement");

        ClearAndResize(_slice_displacement, _slices.size());
        ClearAndResize(_slice_tx, _slices.size());
        ClearAndResize(_slice_ty, _slices.size());
        ClearAndResize(_slice_tz, _slices.size());

        double disp_sum_total = 0;
        int num_voxel_total = 0;

        #pragma omp parallel for reduction(+: disp_sum_total, num_voxel_total)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            int num_voxel_slice = 0;
            double slice_disp = -1, disp_sum_slice = 0;
            double tx_slice = 0, ty_slice = 0, tz_slice = 0;
            double tx_sum_slice = 0, ty_sum_slice = 0, tz_sum_slice = 0;

            if (!_slice_excluded[inputIndex]) {
                // candidate for VoxelFunction
                for (int l = 0; l < _slices[inputIndex].GetT(); l++) {
                    for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
                        for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                            if (InSliceMask(inputIndex,i, j)) {
                                double x = i, y = j, z = 0;
                                _slices[inputIndex].ImageToWorld(x, y, z);
                                double tx = x, ty = y, tz = z;
                                _transformations[inputIndex].Transform(tx, ty, tz);
                                disp_sum_slice += sqrt((tx - x) * (tx - x) + (ty - y) * (ty - y) + (tz - z) * (tz - z));
                                tx_sum_slice += tx - x;
                                ty_sum_slice += ty - y;
                                tz_sum_slice += tz - z;
                                num_voxel_slice++;
                            }
                        }
                    }
                }
                if (num_voxel_slice > 0) {
                    slice_disp = disp_sum_slice / num_voxel_slice;
                    tx_slice = tx_sum_slice / num_voxel_slice;
                    ty_slice = ty_sum_slice / num_voxel_slice;
                    tz_slice = tz_sum_slice / num_voxel_slice;
                    disp_sum_total += disp_sum_slice;
                    num_voxel_total += num_voxel_slice;
                }
            }

            _slice_displacement[inputIndex] = slice_disp;
            _slice_tx[inputIndex] = tx_slice;
            _slice_ty[inputIndex] = ty_slice;
            _slice_tz[inputIndex] = tz_slice;
        }

        return num_voxel_total > 0 ? disp_sum_total / num_voxel_total : -1;
    }

    // -----------------------------------------------------------------------------
    // Calculate Weighted Displacement
    // -----------------------------------------------------------------------------
    double ReconstructionModel::CalculateWeightedDisplacement() {
        logger::debug("CalculateWeightedDisplacement");

        ClearAndResize(_slice_weighted_displacement, _slices.size());

        double disp_sum_total = 0, weight_total = 0;

        #pragma omp parallel for reduction(+: disp_sum_total, weight_total)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            double disp_sum_slice = 0, weight_slice = 0, slice_disp = -1;

            if (!_slice_excluded[inputIndex]) {
                // candidate for VoxelFunction
                for (int l = 0; l < _slices[inputIndex].GetT(); l++) {
                    for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
                            if (InSliceMask(inputIndex,i, j)) {
                                double x = i, y = j, z = 0;
                                _slices[inputIndex].ImageToWorld(x, y, z);
                                double tx = x, ty = y, tz = z;
                                _transformations[inputIndex].Transform(tx, ty, tz);
                                disp_sum_slice += _slice_weight[inputIndex] * _weights[inputIndex](i, j, 0, l) 
                                                    * sqrt((tx - x) * (tx - x) + (ty - y) * (ty - y) + (tz - z) * (tz - z));
                                weight_slice += _slice_weight[inputIndex] * _weights[inputIndex](i, j, 0, l);
                            }
                        }
                    }
                }
                if (weight_slice > 0) {
                    slice_disp = disp_sum_slice / weight_slice;
                    disp_sum_total += disp_sum_slice;
                    weight_total += weight_slice;
                }
            }

            _slice_weighted_displacement[inputIndex] = slice_disp;
        }

        return weight_total > 0 ? disp_sum_total / weight_total : -1;
    }

    // -----------------------------------------------------------------------------
    // Calculate Target Registration Error
    // -----------------------------------------------------------------------------
    double ReconstructionModel::CalculateTRE() {
        logger::debug("CalculateTRE");

        ClearAndReserve(_slice_tre, _slices.size());

        double tre_sum_total = 0;
        int num_voxel_total = 0;

        #pragma omp parallel for reduction(+: tre_sum_total, num_voxel_total)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            double tre_sum_slice = 0, slice_tre = -1;
            int num_voxel_slice = 0;

            if (!_slice_excluded[inputIndex]) {
                // candidate for VoxelFunction
                for (int l = 0; l < _slices[inputIndex].GetT(); l++) {
                    for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
                        for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                            if (InSliceMask(inputIndex,i, j)) {
                                double x = i, y = j, z = 0;
                                _slices[inputIndex].ImageToWorld(x, y, z);
                                double cx = x, cy = y, cz = z;
                                double px = x, py = y, pz = z;
                                _transformations[inputIndex].Transform(cx, cy, cz);
                                _ref_transformations[inputIndex].Transform(px, py, pz);
                                tre_sum_slice += sqrt((cx - px) * (cx - px) + (cy - py) * (cy - py) + (cz - pz) * (cz - pz));
                                num_voxel_slice++;
                            }
                        }
                    }
                }

                if (num_voxel_slice > 0) {
                    slice_tre = tre_sum_slice / num_voxel_slice;
                    tre_sum_total += tre_sum_slice;
                    num_voxel_total += num_voxel_slice;
                }
            }

            _slice_tre[inputIndex] = slice_tre;
        }

        return num_voxel_total > 0 ? tre_sum_total / num_voxel_total : -1;
    }

    // -----------------------------------------------------------------------------
    // Smooth Transformations
    // -----------------------------------------------------------------------------
    void ReconstructionModel::SmoothTransformations(double sigma_seconds, int niter, bool use_slice_inside) {
        logger::debug("SmoothTransformations\n\tsigma = " + to_string(sigma_seconds) + " s");

        //Reset origin for transformations
        GreyImage t = _reconstructed;
        RigidTransformation offset;
        ResetOrigin(t, offset);
        if (_debug)
            offset.Write("reset_origin.dof");

        const Matrix& mo = offset.GetMatrix();
        const Matrix imo = mo.Inverse();
        for (size_t i = 0; i < _transformations.size(); i++)
            _transformations[i].PutMatrix(imo * _transformations[i].GetMatrix() * mo);

        //initial weights
        Matrix parameters(6, _transformations.size());
        Matrix weights(6, _transformations.size());
        ofstream fileOut("motion.txt", ofstream::out | ofstream::app);

        for (size_t i = 0; i < _transformations.size(); i++) {
            parameters(0, i) = _transformations[i].GetTranslationX();
            parameters(1, i) = _transformations[i].GetTranslationY();
            parameters(2, i) = _transformations[i].GetTranslationZ();
            parameters(3, i) = _transformations[i].GetRotationX();
            parameters(4, i) = _transformations[i].GetRotationY();
            parameters(5, i) = _transformations[i].GetRotationZ();

            //write unprocessed parameters to file
            for (int j = 0; j < 6; j++) {
                fileOut << parameters(j, i);
                if (j < 5)
                    fileOut << ",";
                else
                    fileOut << endl;
            }

            for (int j = 0; j < 6; j++)
                weights(j, i) = 0;             //initialise as zero

            if (!_slice_excluded[i]) {
                if (!use_slice_inside) {    //set weights based on image intensities
                    RealPixel smin, smax;
                    _slices[i].GetMinMax(&smin, &smax);
                    if (smax > -1)
                        for (int j = 0; j < 6; j++)
                            weights(j, i) = 1;
                } else {    //set weights based on _slice_inside
                    if (!_slice_inside.empty() && _slice_inside[i])
                        for (int j = 0; j < 6; j++)
                            weights(j, i) = 1;
                }
            }
        }

        //initialise
        Matrix den(6, _transformations.size());
        Matrix num(6, _transformations.size());
        Matrix kr(6, _transformations.size());

        int nloc = 0;
        for (size_t i = 0; i < _transformations.size(); i++)
            nloc = max(nloc, _loc_index[i] + 1);
        logger::debug("\tnumber of slice-locations = " + to_string(nloc));
        const int dim = _transformations.size() / nloc; // assuming equal number of dynamic images for every slice-location

        //step size for sampling volume in error calculation in kernel regression
        constexpr double nstep = 15;
        int step = ceil(_reconstructed.GetX() / nstep);
        if (step > ceil(_reconstructed.GetY() / nstep))
            step = ceil(_reconstructed.GetY() / nstep);
        if (step > ceil(_reconstructed.GetZ() / nstep))
            step = ceil(_reconstructed.GetZ() / nstep);

        //kernel regression
        Array<double> kernel, error, tmp;
        error.reserve(_transformations.size());
        tmp.reserve(_transformations.size());
        for (int iter = 0; iter < niter; iter++) {
            for (int loc = 0; loc < nloc; loc++) {
                //gaussian kernel for current slice-location
                const double sigma = ceil(sigma_seconds / _slice_dt[loc * dim]);
                for (int j = -3 * sigma; j <= 3 * sigma; j++)
                    kernel.push_back(exp(-(j * _slice_dt[loc * dim] / sigma_seconds) * (j * _slice_dt[loc * dim] / sigma_seconds)));

                //kernel-weighted summation
                #pragma omp parallel for collapse(2)
                for (int par = 0; par < 6; par++) {
                    for (int i = loc * dim; i < (loc + 1) * dim; i++) {
                        if (!_slice_excluded[i]) {
                            for (int j = -3 * sigma; j <= 3 * sigma; j++)
                                if (((i + j) >= loc * dim) && ((i + j) < (loc + 1) * dim)) {
                                    num(par, i) += parameters(par, i + j) * kernel[j + 3 * sigma] * weights(par, i + j);
                                    den(par, i) += kernel[j + 3 * sigma] * weights(par, i + j);
                                }
                        } else {
                            num(par, i) = parameters(par, i);
                            den(par, i) = 1;
                        }
                    }
                }

                kernel.clear();
            }

            //kernel-weighted normalisation
            #pragma omp parallel for collapse(2)
            for (int par = 0; par < 6; par++)
                for (size_t i = 0; i < _transformations.size(); i++)
                    kr(par, i) = num(par, i) / den(par, i);

            //recalculate weights using target registration error with original transformations as targets
            for (size_t i = 0; i < _transformations.size(); i++) {
                RigidTransformation processed;
                processed.PutTranslationX(kr(0, i));
                processed.PutTranslationY(kr(1, i));
                processed.PutTranslationZ(kr(2, i));
                processed.PutRotationX(kr(3, i));
                processed.PutRotationY(kr(4, i));
                processed.PutRotationZ(kr(5, i));

                RigidTransformation orig = _transformations[i];

                //need to convert the transformations back to the original coordinate system
                orig.PutMatrix(mo * orig.GetMatrix() * imo);
                processed.PutMatrix(mo * processed.GetMatrix() * imo);

                int n = 0;
                double e = 0;

                if (!_slice_excluded[i]) {
                    #pragma omp parallel for collapse(3) reduction(+: n, e)
                    for (int ii = 0; ii < _reconstructed.GetX(); ii += step)
                        for (int jj = 0; jj < _reconstructed.GetY(); jj += step)
                            for (int kk = 0; kk < _reconstructed.GetZ(); kk += step)
                                if (_reconstructed(ii, jj, kk, 0) > -1) {
                                    double x = ii, y = jj, z = kk;
                                    _reconstructed.ImageToWorld(x, y, z);
                                    double xx = x, yy = y, zz = z;
                                    orig.Transform(x, y, z);
                                    processed.Transform(xx, yy, zz);
                                    x -= xx;
                                    y -= yy;
                                    z -= zz;
                                    e += sqrt(x * x + y * y + z * z);
                                    n++;
                                }
                }

                if (n > 0) {
                    e /= n;
                    error.push_back(e);
                    tmp.push_back(e);
                } else
                    error.push_back(-1);
            }

            sort(tmp.begin(), tmp.end());
            const double median = tmp[round(tmp.size() * 0.5) - 1];

            if (_debug) {
                stringstream ss;
                if (iter == 0)
                    ss << "\titeration:median_error(mm)...";
                ss << iter << ":" << median << ", ";
                logger::debug(ss.str());
            }

            #pragma omp parallel for
            for (size_t i = 0; i < _transformations.size(); i++) {
                double value = 0;
                if (error[i] >= 0 && !_slice_excluded[i]) {
                    if (error[i] <= median * 1.35)
                        value = 1;
                    else
                        value = median * 1.35 / error[i];
                }
                for (int par = 0; par < 6; par++)
                    weights(par, i) = value;
            }

            error.clear();
            tmp.clear();
        }

        ofstream fileOut2("motion-processed.txt", ofstream::out | ofstream::app);
        ofstream fileOut3("weights.txt", ofstream::out | ofstream::app);
        ofstream fileOut4("outliers.txt", ofstream::out | ofstream::app);
        ofstream fileOut5("empty.txt", ofstream::out | ofstream::app);

        for (size_t i = 0; i < _transformations.size(); i++) {
            fileOut3 << weights(0, i) << " ";

            if (weights(0, i) <= 0)
                fileOut5 << i << " ";

            _transformations[i].PutTranslationX(kr(0, i));
            _transformations[i].PutTranslationY(kr(1, i));
            _transformations[i].PutTranslationZ(kr(2, i));
            _transformations[i].PutRotationX(kr(3, i));
            _transformations[i].PutRotationY(kr(4, i));
            _transformations[i].PutRotationZ(kr(5, i));

            for (int j = 0; j < 6; j++) {
                fileOut2 << kr(j, i);
                if (j < 5)
                    fileOut2 << ",";
                else
                    fileOut2 << endl;
            }
        }

        //Put origin back
        #pragma omp parallel for
        for (size_t i = 0; i < _transformations.size(); i++)
            _transformations[i].PutMatrix(mo * _transformations[i].GetMatrix() * imo);
    }

    // -----------------------------------------------------------------------------
    // Scale Transformations
    // -----------------------------------------------------------------------------
    void ReconstructionModel::ScaleTransformations(double scale) {
        logger::trace((boost::format("Scaling transformations. scale = %1%.") % scale).str());

        if (scale == 1)
            return;

        if (scale < 0)
            throw runtime_error("Scaling of transformations undefined for scale < 0.");

        //Reset origin for transformations
        GreyImage t = _reconstructed;
        RigidTransformation offset;
        ResetOrigin(t, offset);
        if (_debug)
            offset.Write("reset_origin.dof");

        const Matrix& mo = offset.GetMatrix();
        const Matrix imo = mo.Inverse();

        #pragma omp parallel for
        for (size_t i = 0; i < _transformations.size(); i++) {
            _transformations[i].PutMatrix(imo * _transformations[i].GetMatrix() * mo);

            //Scale transformations
            const Matrix orig = logm(_transformations[i].GetMatrix());
            Matrix scaled = orig;
            for (int row = 0; row < 4; row++)
                for (int col = 0; col < 4; col++)
                    scaled(row, col) = scale * orig(row, col);
            _transformations[i].PutMatrix(expm(scaled));

            //Put origin back
            _transformations[i].PutMatrix(mo * _transformations[i].GetMatrix() * imo);
        }
    }

    // -----------------------------------------------------------------------------
    // Super-Resolution of Volume
    // -----------------------------------------------------------------------------
    void ReconstructionModel::Superresolution(int iter) {
        START_TIMING();

        logger::trace((boost::format("Superresolution %1%") % iter).str());

        //Remember current reconstruction for edge-preserving smoothing
        RealImage original = _reconstructed;
        Parallel::Superresolution parallelSuperresolution(this);
        parallelSuperresolution();

        RealImage& addon = parallelSuperresolution.addon;
        _confidence_map = move(parallelSuperresolution.confidence_map);

        if (_verbose) {
            SaveReconImageFile(_confidence_map,"confidence-map",iter);
            SaveReconImageFile(addon,"addon",iter);
        }

        if (!_adaptive) {
            RealPixel *pa = addon.Data();
            RealPixel *pcm = _confidence_map.Data();
            #pragma omp parallel for
            for (int i = 0; i < addon.NumberOfVoxels(); i++) {
                if (pcm[i] > 0) {
                    // ISSUES if pcm[i] is too small leading to bright pixels
                    pa[i] /= pcm[i];
                    //this is to revert to normal (non-adaptive) regularisation
                    pcm[i] = 1;
                }
                //else {
                    //pcm[i] = 1;
                    // Future potential improvement - what should happen in this scenario?
                //}
            }
        }

        _reconstructed += addon * _alpha; //_average_volume_weight;

        //bound the intensities
        RealPixel *pr = _reconstructed.Data();
        RealPixel *pm = _confidence_map.Data();
        //RealPixel *pm = _mask.Data(); //Future potential improvement
        #pragma omp parallel for
        for (int i = 0; i < _reconstructed.NumberOfVoxels(); i++) {
            if (pm[i]>0)
            {
                if (pr[i] < _min_intensity * 0.9)
                    pr[i] = _min_intensity * 0.9;
                if (pr[i] > _max_intensity * 1.1)
                    pr[i] = _max_intensity * 1.1;
            }
            else
                pr[i]=0;
        }
        //Smooth the reconstructed image
        AdaptiveRegularization(iter, original);

        END_TIMING("Superresolution");
    }

    // -----------------------------------------------------------------------------
    // Adaptive Regularization
    // -----------------------------------------------------------------------------
    void ReconstructionModel::AdaptiveRegularization(int iter, const RealImage& original) {
        logger::trace("AdaptiveRegularization.");

        Array<double> factor(13);
        for (int i = 0; i < 13; i++) {
            for (int j = 0; j < 3; j++)
                factor[i] += fabs(double(_directions[i][j]));
            factor[i] = 1 / factor[i];
        }

        Array<RealImage> b(13, _reconstructed);

        Parallel::AdaptiveRegularizationStepA arStepA(this, b, factor, original);
        ParallelForEachVoxel(_reconstructed.Attributes(), _reconstructed, arStepA);

        RealImage original2 = _reconstructed;
        Parallel::AdaptiveRegularizationStepB arStepB(this, b, original2);
        ParallelForEachVoxel(_reconstructed.Attributes(), _reconstructed, arStepB);

        if (_alpha * _lambda / (_delta * _delta) > 0.068)
            cerr << "Warning: regularization might not have smoothing effect! Ensure that alpha*lambda/delta^2 is below 0.068." << endl;
    }

    // -----------------------------------------------------------------------------
    // ReadRefTransformations
    // -----------------------------------------------------------------------------
    void ReconstructionModel::ReadRefTransformations(const char *folder) {
        stringstream ss;
        ss << "Reading reference transformations from: " << folder;
        logger::info(ss.str());
        ReadTransformations(folder, _slices.size(), _ref_transformations);
    }

    // -----------------------------------------------------------------------------
    // Calculate Entropy
    // -----------------------------------------------------------------------------
    double ReconstructionModel::CalculateEntropy() {
        double sum_x_sq = 0;

        #pragma omp parallel for reduction(+: sum_x_sq)
        for (int i = 0; i < _reconstructed.GetX(); i++)
            for (int j = 0; j < _reconstructed.GetY(); j++)
                for (int k = 0; k < _reconstructed.GetZ(); k++)
                    if (_reconstructed_mask(i, j, k) == 1)
                        for (int f = 0; f < _reconstructed.GetT(); f++) {
                            const double x = _reconstructed(i, j, k, f);
                            sum_x_sq += x * x;
                        }

        const double x_max = sqrt(sum_x_sq);
        double entropy = 0;

        #pragma omp parallel for reduction(+: entropy)
        for (int i = 0; i < _reconstructed.GetX(); i++)
            for (int j = 0; j < _reconstructed.GetY(); j++)
                for (int k = 0; k < _reconstructed.GetZ(); k++)
                    if (_reconstructed_mask(i, j, k) == 1)
                        for (int f = 0; f < _reconstructed.GetT(); f++) {
                            const double x = _reconstructed(i, j, k, f);
                            if (x > 0)
                                entropy += x / x_max * log(x / x_max);
                        }

        return -entropy;
    }

    // -----------------------------------------------------------------------------
    // Common save functions
    // -----------------------------------------------------------------------------

    void ReconstructionModel::SaveAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, const string& name, int iter, int rec_iter) {
        logger::trace("Save "+name);

        Array<RealImage> save_stacks;
        save_stacks.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++)
            save_stacks.push_back(RealImage(stacks[i].Attributes()));

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < source.size(); inputIndex++)
            for (int t = 0; t < source[inputIndex].GetT(); t++)
                for (int j = 0; j < source[inputIndex].GetY(); j++)
                    for (int i = 0; i < source[inputIndex].GetX(); i++)
                        save_stacks[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], t) = source[inputIndex](i, j, 0, t);


        #pragma omp parallel for
        for (size_t i = 0; i < stacks.size(); i++) {
            string filename;
                SaveStackFile(save_stacks[i],name,i,iter,rec_iter);
        }
    }

    void ReconstructionModel::Save2dAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, const string& name, int iter, int rec_iter) {
        logger::trace("Save "+name);

        Array<RealImage> save_stacks;
        save_stacks.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++) {
            ImageAttributes attr = stacks[i].Attributes();
            attr._t = 1; //enforce no time.
            save_stacks.push_back(RealImage(attr));
        }

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < source.size(); inputIndex++)
            for (int i = 0; i < source[inputIndex].GetX(); i++)
                for (int j = 0; j < source[inputIndex].GetY(); j++)
                        save_stacks[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], 0) = source[inputIndex](i, j, 0);

        #pragma omp parallel for
        for (size_t i = 0; i < stacks.size(); i++) {
            string filename;
            // if (iter < 0 || rec_iter < 0)
            //     SaveStackFile(save_stacks[i],name,i);
            // else
                SaveStackFile(save_stacks[i],name,i,iter,rec_iter);
        }
    }


    void ReconstructionModel::GetAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, Array<RealImage>& output){
        logger::trace("GetAsStack");

        output.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++)
            output.push_back(RealImage(stacks[i].Attributes()));

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < source.size(); inputIndex++)
            for (int t = 0; t < source[inputIndex].GetT(); t++)
                for (int j = 0; j < source[inputIndex].GetY(); j++)
                    for (int i = 0; i < source[inputIndex].GetX(); i++)
                        output[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], t) = source[inputIndex](i, j, 0, t);
    }
    
    void ReconstructionModel::Get2dAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, Array<RealImage>& output){
        logger::trace("GetAsStack");

        output.reserve(stacks.size());

        for (size_t i = 0; i < stacks.size(); i++)
            output.push_back(RealImage(stacks[i].Attributes()));

        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < source.size(); inputIndex++)
            for (int j = 0; j < source[inputIndex].GetY(); j++)
                for (int i = 0; i < source[inputIndex].GetX(); i++)
                    output[_stack_index[inputIndex]](i, j, _stack_loc_index[inputIndex], 0) = source[inputIndex](i, j, 0, 0);

    }

    /// Save corrected slices
    void ReconstructionModel::SaveCorrectedSlices() {
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _corrected_slices.size(); inputIndex++)
            SaveSliceFile(_corrected_slices[inputIndex], "corrected", inputIndex);
    }

    void ReconstructionModel::SaveCorrectedSlicesAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        SaveAsStack(_corrected_slices, stacks, "correctedstack", iter, rec_iter);
    }

    /// Save error
    void ReconstructionModel::SaveError() {
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _error.size(); inputIndex++)
            SaveSliceFile(_error[inputIndex], "error", inputIndex);
    }


    void ReconstructionModel::SaveErrorAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        SaveAsStack(_error, stacks, "errorstack", iter, rec_iter);
    }
    

    void ReconstructionModel::SaveBiasFields() {
        #pragma omp parallel for
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            SaveSliceFile(_bias[inputIndex], "bias", inputIndex);
    }

    /// Save bias fields
    void ReconstructionModel::SaveBiasFieldsAsStack(const Array<RealImage>& stacks, int iter, int rec_iter) {
        Save2dAsStack(_bias, stacks, "bias", iter, rec_iter);
    }


    void ReconstructionModel::GetBiasFieldsAsStack(const Array<RealImage>& stacks, vector<RealImage>& output){
        Get2dAsStack(_bias, stacks, output);
    }

    // -----------------------------------------------------------------------------
    // Save slice info
    // -----------------------------------------------------------------------------
    void ReconstructionModel::SaveSliceInfo(const char *filename, const Array<string>& stack_files) {
        ofstream info(filename);

        info << setprecision(3);

        // header
        info << "StackIndex" << "\t"
            << "StackLocIndex" << "\t"
            << "LocIndex" << "\t"
            << "InputIndex" << "\t"
            << "File" << "\t"
            << "Scale" << "\t"
            << "StackFactor" << "\t"
            << "Time" << "\t"
            << "TemporalResolution" << "\t"
            << "Included" << "\t" // Included slices
            << "Excluded" << "\t"  // Excluded slices
            << "Outside" << "\t"  // Outside slices
            << "Weight" << "\t"
            << "MeanDisplacement" << "\t"
            << "MeanDisplacementX" << "\t"
            << "MeanDisplacementY" << "\t"
            << "MeanDisplacementZ" << "\t"
            << "WeightedMeanDisplacement" << "\t"
            << "TRE" << "\t"
            << "TranslationX" << "\t"
            << "TranslationY" << "\t"
            << "TranslationZ" << "\t"
            << "RotationX" << "\t"
            << "RotationY" << "\t"
            << "RotationZ" << "\t";
        info << "\b" << endl;

        for (size_t i = 0; i < _slices.size(); i++) {
            RigidTransformation& t = _transformations[i];
            info << (_stack_index.size() > i ? _stack_index[i] : -1 ) << "\t"
                << (_stack_loc_index.size() > i ? _stack_loc_index[i] : -1 ) << "\t"
                << (_loc_index.size() > i ? _loc_index[i] : -1 )  << "\t"
                << i << "\t"
                << (_stack_index.size() > i && stack_files.size() > i ? stack_files[_stack_index[i]] : "" )  << "\t"
                << (_scale.size() > i ? _scale[i] : -1 )  << "\t"
                << (_stack_index.size() > i && _stack_factor.size() > i ? _stack_factor[_stack_index[i]] : -1 )  << "\t"
                << (_slice_dt.size() > i ? _slice_dt[i] : -1 )  << "\t"
                << ((_slice_weight.size() > i && _slice_weight[i] >= 0.5 && _slice_inside.size() > i && _slice_inside[i]) ? 1 : 0) << "\t" // Included slices
                << ((_slice_weight.size() > i && _slice_weight[i] < 0.5 && _slice_inside.size() > i && _slice_inside[i]) ? 1 : 0) << "\t"  // Excluded slices
                << ((_slice_inside.size() > i && !(_slice_inside[i])) ? 1 : 0) << "\t"  // Outside slices
                << (_slice_weight.size() > i ? _slice_weight[i]  : -1 )  << "\t"
                << (_slice_displacement.size() > 0 ? _slice_displacement[i]  : -1 )  << "\t"   
                << (_slice_tx.size() > i ? _slice_tx[i]  : -1 )  << "\t"   
                << (_slice_ty.size() > i ? _slice_ty[i]  : -1 )  << "\t"   
                << (_slice_tz.size() > i ? _slice_tz[i]  : -1 )  << "\t"   
                << (_slice_weighted_displacement.size() > i ? _slice_weighted_displacement[i]  : -1 )  << "\t"
                << (_slice_tre.size() > i ? _slice_tre[i]  : -1 )  << "\t"
                << t.GetTranslationX() << "\t"
                << t.GetTranslationY() << "\t"
                << t.GetTranslationZ() << "\t"
                << t.GetRotationX() << "\t"
                << t.GetRotationY() << "\t"
                << t.GetRotationZ() << "\t";
            info << "\b" << endl;
        }
    }

    // -----------------------------------------------------------------------------


    // scale the reconstructed volume
    void ReconstructionModel::ScaleVolume(RealImage& reconstructed) {
        START_TIMING();
        double scalenum = 0, scaleden = 0;
        double sum_weight, sum_slice, sum_sim, sum_sim_weight, sum_slice_weight = 0;

        //#pragma omp parallel for 
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            // alias for the current slice
            const RealImage& slice = _slices[inputIndex];
            const RealImage& mask = _slice_masks[inputIndex];

            //alias for the current weight image
            const RealImage& w = _weights[inputIndex];

            // alias for the current simulated slice
            const RealImage& sim = _simulated_slices[inputIndex];

            for (int l = 0; l < slice.GetT(); l++)
                for (int j = 0; j < slice.GetY(); j++)
                    for (int i = 0; i < slice.GetX(); i++)
                        if (InSliceMask(inputIndex,i, j)) {
                            //scale - intensity matching
                            if (_simulated_weights[inputIndex](i, j, 0, l) > 0.99) {
                                scalenum += w(i, j, 0, l) * _slice_weight[inputIndex] * slice(i, j, 0, l) * sim(i, j, 0, l);
                                scaleden += w(i, j, 0, l) * _slice_weight[inputIndex] * sim(i, j, 0, l) * sim(i, j, 0, l);
                                sum_weight += w(i, j, 0, l);
                                sum_slice += slice(i, j, 0, l);
                                sum_sim += sim(i, j, 0, l);
                                sum_slice_weight += _slice_weight[inputIndex];
                                sum_sim_weight += _simulated_weights[inputIndex](i, j, 0, l);
                            }
                        }

        } //end of loop for a slice inputIndex

        //calculate scale for the volume
        double scale = 1;

        if (scaleden > 0)
            scale = scalenum / scaleden;

        logger::trace((boost::format("scale : %1%") % scale).str());
        logger::trace((boost::format("scalenum : %1%") % scalenum).str());
        logger::trace((boost::format("scaleden : %1%") % scaleden).str());
        logger::trace((boost::format("sum_weight : %1%") % sum_weight).str());
        logger::trace((boost::format("sum_slice : %1%") % sum_slice).str());
        logger::trace((boost::format("sum_sim : %1%") % sum_sim).str());
        logger::trace((boost::format("sum_sim_weight : %1%") % sum_sim_weight).str());
        logger::trace((boost::format("sum_slice_weight : %1%") % sum_slice_weight).str());

        RealPixel *ptr = reconstructed.Data();
        #pragma omp parallel for
        for (int i = 0; i < reconstructed.NumberOfVoxels(); i++)
            if (ptr[i] > 0)
                ptr[i] *= scale;
        END_TIMING("ScaleVolume");
    }


    // initialise slice EM step
    void ReconstructionModel::InitializeEM() {
        START_TIMING();

        ClearAndReserve(_weights, _slices.size());
        ClearAndReserve(_bias, _slices.size());
        ClearAndReserve(_scale, _slices.size());
        ClearAndReserve(_slice_weight, _slices.size());

        for (size_t i = 0; i < _slices.size(); i++) {
            //Create images for voxel weights and bias fields
            const ImageAttributes& attr = _slices[i].Attributes();
            _weights.push_back(_slices[i]);
            //_bias.push_back(_slices[i].GetRegion(0,0,0,0,attr._x, attr._y, attr._z, 1)); //Bias 3D only: take first time point 
            _bias.push_back(_slices[i].GetFrame(0,0)); //Bias 3D only: take first time point 

            //Create and initialize scales
            _scale.push_back(1);

            //Create and initialize slice weights
            _slice_weight.push_back(1);
        }

        //Find the range of intensities
        _max_intensity = voxel_limits<RealPixel>::min();
        _min_intensity = voxel_limits<RealPixel>::max();

        #pragma omp parallel for reduction(min: _min_intensity) reduction(max: _max_intensity)
        for (size_t i = 0; i < _slices.size(); i++) {
            RealImage& slice = _slices[i];
            for (int t = 0; t < slice.GetT(); ++t){
                for (int y = 0; y < slice.GetY(); ++y){
                    for (int x = 0; x < slice.GetX(); ++x) {
                        if (InSliceMask(i,x,y) && slice(x,y,0,t) > 0){
                            if (slice(x,y,0,t) > _max_intensity)
                                _max_intensity = slice(x,y,0,t);
                            if (slice(x,y,0,t) < _min_intensity)
                                _min_intensity = slice(x,y,0,t);
                        } 
                    }
                }
            }
        }
        END_TIMING("InitializeEM");
    }



    // initialise parameters of EM robust statistics
    void ReconstructionModel::InitializeRobustStatistics() {
        START_TIMING();
        Array<int> sigma_numbers(_slices.size());
        Array<double> sigma_values(_slices.size());
        RealImage slice;

        #pragma omp parallel for private(slice)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            slice = _slices[inputIndex];
            const RealImage& mask = _slice_masks[inputIndex];
            //Voxel-wise sigma will be set to stdev of volumetric errors
            //For each slice voxel
            for (int l = 0; l < slice.GetT(); l++)
                for (int i = 0; i < slice.GetX(); i++)
                    for (int j = 0; j < slice.GetY(); j++)
                        if (InSliceMask(inputIndex,i,j)) { 
                            //calculate stev of the errors
                            if (_simulated_inside[inputIndex](i, j, 0, l) == 1 && _simulated_weights[inputIndex](i, j, 0, l) > 0.99) {
                                slice(i, j, 0, l) -= _simulated_slices[inputIndex](i, j, 0, l);
                                sigma_values[inputIndex] += slice(i, j, 0, l) * slice(i, j, 0, l);
                                sigma_numbers[inputIndex]++;
                            }
                        }

            //if slice does not have an overlap with ROI, set its weight to zero
            if (!_slice_inside[inputIndex])
                _slice_weight[inputIndex] = 0;
        }

        double sigma = 0;
        int num = 0;
        #pragma omp simd
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            sigma += sigma_values[inputIndex];
            num += sigma_numbers[inputIndex];
        }

        //Force exclusion of slices predefined by user
        for (size_t i = 0; i < _force_excluded.size(); i++)
            if (_force_excluded[i] > 0 && _force_excluded[i] < _slices.size())
                _slice_weight[_force_excluded[i]] = 0;

        //initialize sigma for voxel-wise robust statistics
        _sigma = sigma / num;

        //initialize sigma for slice-wise robust statistics
        _sigma_s = 0.025;
        //initialize mixing proportion for inlier class in voxel-wise robust statistics
        _mix = 0.9;
        //initialize mixing proportion for outlier class in slice-wise robust statistics
        _mix_s = 0.9;
        //Initialise value for uniform distribution according to the range of intensities
        _m = 1 / (2.1 * _max_intensity - 1.9 * _min_intensity);

        logger::trace((boost::format("Initializing robust statistics: sigma=%1% m=%2% mix=%3% mix_s=%4%") % sqrt(_sigma) % _m % _mix % _mix_s).str());
        END_TIMING("InitializeRobustStatistics");
    }

    
    // structure-based outlier rejection step (NCC-based)
    void ReconstructionModel::GlobalStructuralExclusion() {
        START_TIMING();

        double source_padding = -1;
        const double target_padding = -inf;
        constexpr bool dofin_invert = false;
        constexpr bool twod = false;

        RealImage& source = _reconstructed;
        source *= _evaluation_mask;
        RealPixel smin, smax;
        source.GetMinMax(&smin, &smax);

        if (smin < -0.1)
            source_padding = -1;
        else if (smin < 0.1)
            source_padding = 0;

        Array<double> reg_ncc(_slices.size());
        double mean_ncc = 0;
        int number_of_excluded = 0;

        stringstream ss;
        ss << " - excluded : ";

        // #pragma omp parallel for private(output, target, slice_mask) firstprivate(imagetransformation) reduction(+: mean_ncc)
        for (size_t inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
            // transfrom reconstructed volume to the slice space

            RealImage output, target, slice_mask;
            output.Initialize(_slices[inputIndex].Attributes());

            GenericLinearInterpolateImageFunction<RealImage> interpolator;
            ImageTransformation imagetransformation;
            imagetransformation.Input(&source);
            imagetransformation.Output(&output);
            imagetransformation.TargetPaddingValue(target_padding);
            imagetransformation.SourcePaddingValue(source_padding);
            imagetransformation.Interpolator(&interpolator);
            imagetransformation.TwoD(twod);
            imagetransformation.Invert(dofin_invert);
            imagetransformation.Transformation(&_transformations[inputIndex]);

            imagetransformation.Run();
            
            // blur the original slice
            GaussianBlurringWithPadding<RealPixel> gb(target.GetXSize() * 0.6, source_padding);
            RealImage maskedSlice = _slices[inputIndex] * _slice_masks[inputIndex];
            gb.Input(&maskedSlice);
            gb.Output(&target);
            gb.Run();
            
            // compute NCC
            double output_ncc = ComputeNCC(target, output, 0);
            if (output_ncc == -1)
                output_ncc = 1;
            reg_ncc[inputIndex] = output_ncc;
            mean_ncc += output_ncc;

            double current_global_NCC_threshold = _global_NCC_threshold;
            if (_current_iteration == 0)
                current_global_NCC_threshold = _global_NCC_threshold * 0.75;

            // set slice weight
            if (output_ncc > current_global_NCC_threshold) {
                _structural_slice_weight[inputIndex] = 1;
            } else {
                _structural_slice_weight[inputIndex] = -1;
                number_of_excluded++;
                if (_debug)
                    ss << inputIndex << "(" << output_ncc << "), ";
            }

        }
        ss << endl;
        mean_ncc /= _slices.size();

        ss << " - mean registration ncc: " << mean_ncc << " | " << number_of_excluded << " / " << _slices.size() << endl;
        logger::debug(ss.str());

        END_TIMING("GlobalStructuralExclusion");
    }


} // namespace svr

// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR
#include "svr/ReconstructionAlgorithm.h"
#include "svr/Profiling.h"
#include "svr/Parallel.h"

using namespace std;
using namespace mirtk;
using namespace svr::Utility;

namespace svr {

    ReconstructionAlgorithm::ReconstructionAlgorithm(ReconstructionParameters &parameters) 
        : _parameters(parameters) {

        // Verify parameters are valid before creating object.
        parameters.Validate();

        // Set local debug parameters (e.g. profile class looks for these)
        _verbose = _parameters.verbose;
        _debug = _parameters.debug || _parameters.verbose; // verbose makes debug implicit.
        _profile = _parameters.profile;

        ConfigureLogging(_parameters.logID);

        // Write parameters to log
        stringstream ss;
        ss << "Reconstruction 4D input parameters: " << endl;
        ss << _parameters;
        logger::debug(ss.str());
    }

    int ReconstructionAlgorithm::Run() {
        START_TIMING();
        boost::timer::cpu_timer reconstruction_timer; 
        
        if(Initialise())
            return 1;

        reconstruction_timer.start();

        if(PreProcess()) 
            return 1;

        if(Process())  
            return 1;

        reconstruction_timer.stop();

        if (_parameters.iterations > 0) // Only evaluate if more than one iteration executed.
            Evaluate(reconstruction_timer);

        if(PostProcess()) 
            return 1;

        logger::info("Completed successfully.");

        END_TIMING("Reconstruct::Run");
        return 0;
    }

    int ReconstructionAlgorithm::Initialise() {
        START_TIMING();
        logger::info("Initialising... ");
        logger::flush();

        //Set debuging modes
        ConfigureDebugMode();

        // Initialisation of MIRTK image reader library
        InitializeIOLibrary();

        // Read stack files
        ReadStackFiles();

        //Read stack transformations
        ReadStackTransformationFiles();

        // Number of packages for each stack
        ReadPackagesParam();

        // Set structural flag
        _reconstruction.SetStructural(_parameters.structural);

        END_TIMING("Reconstruction4DFunction::Initialise");
        return 0;
    }

    int ReconstructionAlgorithm::PreProcess()
    {
        START_TIMING();
        logger::info("Pre-processing... ");
        logger::flush();

        CalculateSliceThickness();                                      // Calculate slice thickness per stack

        if (_parameters.packages.size() > 0)
            SplitStacksIntoPackages();                                  // Split stacks into packages where applicable

        _reconstruction.SetForceExcludedSlices(_parameters.forceExcludedSlices); //Set force excluded slices   
        _reconstruction.SetForceExcludedStacks(_parameters.forceExcludedStacks);

        unique_ptr<RealImage> reconstructionMask;
        RealImage maskCropped;
        SetupVolumeMask(reconstructionMask, maskCropped);               // Binary mask for final volume

        if (_parameters.structural)
            InitialiseVolume(*reconstructionMask);
        else
            InitialiseVolume(maskCropped);

        _reconstruction.SetMask(*reconstructionMask, _parameters.smoothMask);
        if (_debug) ReconstructionModel::SaveReconImageFile(_reconstruction.GetReconstructed(),"initialisedAndMasked"); 

        if (_parameters.stackRegFolder.size() == 0)                     //stack-stack registration if not input
            _reconstruction.StackRegistrations(_stacks, _stackTransformations, _parameters.templateNumber);

        if (_parameters.removeBlackBackground)
            RemoveBlackBackground();                                    // Generate mask from black background
                                                    
        RealImage average = _reconstruction.CreateAverage(_stacks, _stackTransformations);  // Create average
        if (_verbose) ReconstructionModel::SaveStackFile(average, "average");

        RealImage templateStack = _stacks[_parameters.templateNumber];  // Copy of uncropped template (for performSecondStackReg)

        CropStacks();

        if (_parameters.performSecondStackReg)                          // Optional second stack-stack registration
           _reconstruction.StackRegistrationsWithTemplate(_stacks, _stackTransformations, templateStack);
           
        if (_parameters.intensityMatchingMode == RESCALE_STACKS_ONLY || 
            _parameters.intensityMatchingMode == FULL_RESCALING)        // Rescale stack intensities
            _reconstruction.MatchStackIntensitiesWithMasking(_stacks, _stackTransformations, _parameters.averageValue);
         
        average = _reconstruction.CreateAverage(_stacks, _stackTransformations);
        if (_verbose) 
            ReconstructionModel::SaveStackFile(average, "average2");    // Second average

        _reconstruction.CreateSlicesAndTransformations(_stacks, _stackTransformations, _parameters.thickness);  //Create slices and transformations
        
        if (_debug)_reconstruction.InitCorrectedSlices();
        if (_debug) _reconstruction.InitError();
        
        if(_parameters.sliceMaskingMode > MASK_MODE::NO_SLICE_MASKS)
            _reconstruction.CreateSliceMasks();                         // Create slice masks

        if(_debug || _parameters.saveSlicesFlag) _reconstruction.SaveSlices();
        if (_debug) _reconstruction.SaveSliceMasks();
        
        if (!_parameters.frameRegfolder.empty())                        //Read slice transformations if provided
            _reconstruction.ReadTransformations(_parameters.frameRegfolder.c_str());        
        else if (!_parameters.sliceLocRegFolder.empty())
            _reconstruction.ReadSliceTransformations(_parameters.sliceLocRegFolder.c_str());  

        _reconstruction.InitTRE();                                      // Initialise TRE

        if (_parameters.sigma > 0) 
            _reconstruction.SetSigma(_parameters.sigma);                //Set sigma for the bias field smoothing
        else _reconstruction.SetSigma(20); 

        _reconstruction.InitializeEM();                                 //Initialise _parameters structures for EM

        END_TIMING("Reconstruction4DFunction::PreProcess");
        return 0;
    }

    int ReconstructionAlgorithm::Process() {
        START_TIMING();
        logger::info("Processing reconstruction...");
        logger::flush();

        Array<double> entropyList;

        for (int iter = 0; iter < _parameters.iterations; iter++) { // REGISTRATION ITERATIONS
            logger::info("------------------------------------------------------");
            logger::info((boost::format("Registration iteration %1% / %2%") % (iter+1) % _parameters.iterations).str());
            _reconstruction.SetCurrentIteration(iter);

            if (_verbose) ReconstructionModel::SaveReconImageFile(_reconstruction.GetReconstructed(),"newIteration",iter);

            if (_parameters.structural)
                _reconstruction.StaticMaskReconstructedVolume(); 
                
            if (_verbose) ReconstructionModel::SaveReconImageFile(_reconstruction.GetReconstructed(),"maskedVolume",iter);

            SliceToVolumeRegistration(iter);                        //perform slice-to-volume registrations

            if (_parameters.structural)
                _reconstruction.GlobalStructuralExclusion();         // Run global NNC structure-based outlier rejection of slices

            SetSmoothingParameters(iter);                           //Set smoothing parameters

            if (iter < _parameters.iterations - 1) 
                _reconstruction.SpeedupOn(); 
            else _reconstruction.SpeedupOff();                      //Use faster _reconstruction during iterations and slower for final _reconstruction

            if (_parameters.robustSlicesOnly) 
                _reconstruction.ExcludeWholeSlicesOnly();           //Exclude whole slices only

            _reconstruction.InitializeEMValues();                   //Initialise values of weights, scales and bias fields

            _reconstruction.CoeffInit();                   //Calculate matrix of transformation between voxels of slices and volume
            if (_verbose) ReconstructionModel::SaveReconImageFile(_reconstruction.GetVolumeWeights(), "weights", iter);

            _reconstruction.GaussianReconstruction();      //Initialize reconstructed image with Gaussian weighted _reconstruction
            if (_verbose) ReconstructionModel::SaveReconImageFile(_reconstruction.GetReconstructed(),"gaussianRecon",iter);
            
            entropyList.clear();                                    // Calculate entropy
            entropyList.push_back(_reconstruction.CalculateEntropy());

            
            _reconstruction.SimulateSlices();              //Simulate slices (needs to be done after Gaussian _reconstruction)
            SaveIntermediateSimulatedSlices(iter);                  //Save intermediate simulated slices
            _reconstruction.InitializeRobustStatistics();           //Initialize robust statistics parameters

            if (_parameters.robustStatistics) _reconstruction.EStep(); //EStep
            if (_verbose) _reconstruction.SaveWeightsAsStack(_stacks, iter);

    
            if (_parameters.structural) {                   // Run local SSIM structure-based outlier rejection
                _reconstruction.SStep();
                _reconstruction.SaveSSIMMapsAsStack(_stacks,iter);

            }
            
            // SR RECONSTRUCTION ITERATIONS
            int srIterations = CalcReconstuctionIterations(iter);                      //number of _reconstruction iterations
            for (int i = 0; i < srIterations; i++) {
                logger::debug("------------------------------------------------------");
                logger::info((boost::format("SR-Reconstruction iteration %1% / %2%") % (i+1) % srIterations).str());

                if (_parameters.intensityMatchingMode >= IM_MODE::FULL_RESCALE_EXCEPT_STACKS) {
                    if (_parameters.sigma > 0) 
                        _reconstruction.Bias();                     //calculate bias fields    
                    _reconstruction.Scale();                        //calculate scales
                }

                _reconstruction.Superresolution(i);        //Update reconstructed volume via superresolution

                if (_parameters.intensityMatchingMode >= IM_MODE::FULL_RESCALE_EXCEPT_STACKS && _parameters.sigma > 0)
                    _reconstruction.NormaliseBias(iter, i);

                
                if (_debug) {                            //Save intermediate reconstructed volume
                    RealImage reconstructed = _reconstruction.GetReconstructed();
                    StaticMaskVolume(reconstructed, _reconstruction.GetMask(), -1);
                    ReconstructionModel::SaveReconImageFile(reconstructed, iter, i );
                }

                entropyList.push_back(
                    _reconstruction.CalculateEntropy());            // Calculate Entropy
                _reconstruction.SimulateSlices();          // Simulate slices (needs to be done after the update of the reconstructed volume)
                
                // If not last iteration
                if (_parameters.lastIterRSCalcs || (!_parameters.lastIterRSCalcs && (i + 1 < _parameters.srIterations))) {

                    SaveIntermediateSimulatedSlices(iter, i);

                    if (_parameters.robustStatistics) {
                        _reconstruction.MStep(i + 1);
                        _reconstruction.EStep();
                    }


                    // Run local SSIM structure-based outlier rejection
                    if (_parameters.structural) {
                        _reconstruction.SStep();
                        _reconstruction.SaveSSIMMapsAsStack(_stacks,iter, i);
                    }

                }

            }// END RECONSTRUCTION ITERATIONS

            _reconstruction.StaticMaskReconstructedVolume();      //Mask reconstructed image to ROI given by the reconstructionMask

            RealImage reconstructed = _reconstruction.GetReconstructed();
            StaticMaskVolume(reconstructed, _reconstruction.GetMask(), -1);  
            if (_debug) ReconstructionModel::SaveReconImageFile(reconstructed, iter);    //Save reconstructed image
            _entropy.push_back(entropyList);                                //Save Calculated Entropy

            CalculateDisplacements(iter, reconstructed);

            if (_verbose) {
                // Display Displacements and TRE
                logger::debug("Mean Displacement (iter "+to_string(iter)+") = "+to_string( _meanDisplacement[iter])+"mm");
                logger::debug("Mean Weighted Displacement (iter "+to_string(iter)+") = "+to_string(_meanWeightedDisplacement[iter])+" mm");

                // Save Info for Iteration
                logger::debug("SaveSliceInfo");
                _reconstruction.SaveSliceInfo((boost::format("info_mc%02i.tsv") % iter).str().c_str(), _parameters.stackFiles);
            }
        }// END REGISTRATION ITERATIONS

        END_TIMING("Reconstruction4DFunction::Process");
        return 0;
    }

    int ReconstructionAlgorithm::PostProcess()
    {
        START_TIMING();
        logger::info("Post-processing... ");
        logger::flush();

        //save final result
        logger::debug("RestoreSliceIntensities");
        if (_parameters.intensityMatchingMode >= IM_MODE::FULL_RESCALE_EXCEPT_STACKS)
            _reconstruction.RestoreSliceIntensities();
        
        logger::debug("ScaleVolume");
        if (_parameters.intensityMatchingMode >= IM_MODE::FULL_RESCALE_EXCEPT_STACKS)
           _reconstruction.ScaleVolume();
        
        logger::debug("Saving Reconstructed Volume");
        _reconstruction.GetReconstructed().Write(_parameters.outputName.c_str());

        if(_verbose) {
            logger::trace("SaveSlices");
            _reconstruction.SaveSlicesAsStack(_stacks);

            logger::trace("SaveTransformations");
            _reconstruction.SaveTransformations();
        }

        if(_parameters.saveSlicesFlag) {
            _reconstruction.SaveSlices();
            _reconstruction.SaveTransformations();
        }


        if (!_parameters.infoFilename.empty() && _verbose) {
            logger::trace("SaveSliceInfo4D");
            _reconstruction.SaveSliceInfo(_parameters.infoFilename.c_str(), _parameters.stackFiles);
        }

        if (_verbose) {
            logger::trace("SaveWeights");
            _reconstruction.SaveWeightsAsStack(_stacks);
            if(_parameters.saveSlicesFlag) 
                _reconstruction.SaveWeights();
            logger::trace("SaveBiasFields");
            _reconstruction.SaveBiasFieldsAsStack(_stacks);
            if(_parameters.saveSlicesFlag) 
                _reconstruction.SaveBiasFields();
            
        }

        if (_debug || _parameters.saveEvaluationData){
            logger::trace("SaveSimulatedSlices");
            _reconstruction.SaveSimulatedSlicesAsStack(_stacks);
        }

        if (_debug){
            if(_parameters.saveSlicesFlag) 
                _reconstruction.SaveSimulatedSlices();
        }


        logger::info((boost::format("Output volume :  %1%") % _parameters.outputName).str());

        logger::debug("Reconstruction complete.");
        END_TIMING("Reconstruction4DFunction::PostProcess");
        return 0;
    }

    void ReconstructionAlgorithm::Evaluate(const boost::timer::cpu_timer& reconstruction_timer){
        START_TIMING();
        logger::info("Evaluating reconstruction...");
        logger::flush();

        stringstream ss;
        ss << endl;
        ss << "------------------------------------------------------" << endl;
        ss << "             RECONSTRUCTION EVALUATION                " << endl;
        ss << "------------------------------------------------------" << endl;
        ss << endl;
        
        double ncc = _reconstruction.CalculateNCC(true);
        double nrmse = _reconstruction.CalculateNRMSE(true);
        double rse = _reconstruction.CalculateRSE(true);
        int includedSlices = 0;
        int excludedSlices = 0;
        int outsideSlices = 0;
        _reconstruction.SliceSummary(includedSlices, excludedSlices, outsideSlices);
        
        ss << "Reconstruction quality metrics: " << endl;
        ss << "Similarity measure (Normalised Cross Correlation): " << ncc << endl;
        ss << "Error measure (Normalised Root Mean Square Error): " << nrmse << endl;
        ss << "Slice use summary: " << includedSlices << " included, " <<  excludedSlices << " excluded, " << outsideSlices << " outside ROI." << endl;  
        ss << "Reconstruction time: " << std::stof(reconstruction_timer.format(3,"%w")) << endl;
        ss << endl;

        ss << "------------------------------------------------------" << endl;
        ss << endl;    

        logger::info(ss.str());

        _reconstruction.SaveMetaData(reconstruction_timer);

        END_TIMING("Reconstruct::Evaluate");
    }

    void ReconstructionAlgorithm::ReadPackagesParam()
    {
        if (!_parameters.packages.empty())
        {
            // Specified a package number per stack
            if (_parameters.packages.size() == _stacks.size())
            {
                stringstream ss;
                ss << "Package number : ";
                for (size_t i = 0; i < _stacks.size(); i++)
                {
                    ss << _parameters.packages[i] << " ";
                }
                ss << endl;
                logger::info(ss.str());
            }
            // Single default number of packages for all stacks
            else if (_parameters.packages.size() == 1)
            {
                stringstream ss;
                ss << "Package number : ";
                for (size_t i = 0; i < _stacks.size(); i++)
                {
                    _parameters.packages.push_back(_parameters.packages[0]);
                    ss << _parameters.packages[i] << " ";
                }
                ss << endl;
                logger::info(ss.str());
            }
        }
    }

    void ReconstructionAlgorithm::ReadStackTransformationFiles()
    {
        if (!_parameters.stackRegFolder.empty())
        {
            for (size_t i = 0; i < _stacks.size(); i++)
            {
                logger::info("Reading transformation " + _parameters.stackRegFolder[i]);
                Transformation *transformation = _parameters.stackRegFolder[i] == "id" ? new RigidTransformation : Transformation::New(_parameters.stackRegFolder[i].c_str());
                unique_ptr<RigidTransformation> rigidTransf(dynamic_cast<RigidTransformation *>(transformation));
                _stackTransformations.push_back(*rigidTransf);
            }
            InvertStackTransformations(_stackTransformations);
        }

        //If transformations were not defined by user, initialise them to identity
        if (_stackTransformations.empty())
            _stackTransformations = Array<RigidTransformation>(_stacks.size());
    }
    void ReconstructionAlgorithm::ReadStackFiles()
    {
        int tCount;
        for (int i = 0; i < _parameters.stackFiles.size(); i++)
        {
            logger::debug("Reading stack " + _parameters.stackFiles[i]);
            RealImage stack;

            UniquePtr<ImageReader> image_reader;
            UniquePtr<BaseImage> tmp_image;
            image_reader.reset(ImageReader::TryNew(_parameters.stackFiles[i].c_str()));
            tmp_image.reset(image_reader->Run());
            stack = *tmp_image;
            _stacks.push_back(move(stack));
        }
        logger::debug("Number of stack time points: " + to_string(tCount));

        if (_debug)
            for (int i = 0; i < _parameters.stackFiles.size(); i++)
                ReconstructionModel::SaveStackFile(_stacks[i], "input", i);
    }

    void ReconstructionAlgorithm::ConfigureDebugMode()
    {
        //Set debug flags
        if (_verbose)
        {
            _reconstruction.VerboseOn();
            _reconstruction.DebugOn();
        }
        else if (_debug)
            _reconstruction.DebugOn();
        else
            _reconstruction.DebugOff();

        logger::debug("Debug mode: " + to_string(_debug));
        logger::debug("Verbose debug mode: " + to_string(_verbose));

        //Set profile mode
        if (_parameters.profile)
            _reconstruction.ProfileOn();
        else
            _reconstruction.ProfileOff();
        logger::debug("Profile mode: " + to_string(_parameters.profile));
    }

    void ReconstructionAlgorithm::CropStacks()
    {
        // Mask is transformed to the all stacks and they are cropped
        for (size_t i = 0; i < _stacks.size(); i++)
        {
            // transform the reconstructionMask
            RealImage m = _reconstruction.GetMask();
            TransformMask(_stacks[i], m, _stackTransformations[i]);

            // Crop template stack
            if (_parameters.noStackZCrop)
                CropImageIgnoreZ(_stacks[i], m);
            else
                CropImage(_stacks[i], m);

            if (_debug)
            {
                ReconstructionModel::SaveStackFile(m, "mask", i);
                ReconstructionModel::SaveStackFile(_stacks[i], "cropped", i);
            }
        }
    }

    void ReconstructionAlgorithm::RemoveBlackBackground()
    {
        logger::debug("Creating mask from black background ... ");
        RealImage newmask = CreateMaskFromBlackBackground(&_reconstruction, _stacks, _stackTransformations);
        // reconstructionMask = unique_ptr<RealImage>(new RealImage(newmask));
        if (_verbose)
            ReconstructionModel::SaveReconImageFile(newmask, "maskFromBlackBackground");
        _reconstruction.SetMask(newmask, _parameters.smoothMask, 0.01);
    }

    void ReconstructionAlgorithm::InitialiseVolume(mirtk::RealImage &maskCropped)
    {
        // Initialise reconstructed volume with isotropic resolution
        // if resolution==0 it will be determined from in-plane resolution of the image
        if (_parameters.resolution <= 0)
            _parameters.resolution = _reconstruction.GetReconstructedResolutionFromTemplateStack(_stacks[_parameters.templateNumber]);
        if (_debug)
            logger::info("Initialising volume with isotropic voxel size " + to_string(_parameters.resolution) + "mm");

        // Set temporal resolution
        int temporalSize = _stacks[_parameters.templateNumber].GetT();
        double temporalResolution;
        if (_parameters.timeResolution > 0)
        {
            temporalResolution = _parameters.timeResolution;
        }
        else
        {
            temporalResolution = 0;
            for (auto stack : _stacks)
                temporalResolution += stack.TSize() / _stacks.size();
            logger::debug("Calculated temporalResolution as average of stacks: " + to_string(temporalResolution));
        }


        if (!_parameters.structural) {

        // Original 4D behaviour - create template from mask
        if (_verbose)
            ReconstructionModel::SaveReconImageFile(maskCropped, "templateCropped");
        _reconstruction.CreateTemplateFromStaticMask(maskCropped, _parameters.resolution,
                                                     temporalSize, temporalResolution);

        }
        else {
            // Original 3D behaviour (create template from first stack, required for global structural exclusion)

            RealImage maskedTemplate;
            // Before creating the template we will crop template stack according to the given mask

                // First resample the mask to the space of the stack
            // For template stack the transformation is identity
            RealImage m = maskCropped;
            TransformMask(_stacks[_parameters.templateNumber], m, _stackTransformations[_parameters.templateNumber]);

            // Crop template stack and prepare template for global volumetric registration
            maskedTemplate = _stacks[_parameters.templateNumber] * m;
                        
            CropImage(maskedTemplate, m);

            if (_verbose) {
                ReconstructionModel::SaveReconImageFile(m, "templateMask");
                ReconstructionModel::SaveReconImageFile(maskedTemplate, "templateCropped");
            }

            _reconstruction.CreateTemplate(maskCropped.Attributes(), 
                                        maskedTemplate, 
                                        _parameters.resolution,
                                        temporalSize, 
                                        temporalResolution);
        }

        if (_debug) ReconstructionModel::SaveReconImageFile(_reconstruction.GetReconstructed(),"initialised"); 

    }

    void ReconstructionAlgorithm::SetupVolumeMask(std::unique_ptr<mirtk::RealImage> &reconstructionMask, mirtk::RealImage &maskCropped)
    {
        if (!_parameters.maskFile.empty())
        {
            logger::info("Applying mask " + _parameters.maskFile);
            reconstructionMask = unique_ptr<RealImage>(new RealImage(_parameters.maskFile.c_str()));
        }
        else
        {
            // If no mask was given - try to create mask from the template image in case it was padded
            // Note - parameter verification only allows code to proceed with no mask if removeBlackBackground set.
            reconstructionMask = unique_ptr<RealImage>(new RealImage(_stacks[_parameters.templateNumber]));
            *reconstructionMask = CreateMask(*reconstructionMask);
            logger::error("Warning : no mask was provided");
        }

        if (_verbose) ReconstructionModel::SaveReconImageFile(*reconstructionMask,"rawmask"); 

        // Initialise Reconstructed Volume
        // Crop reconstructionMask
        maskCropped = *reconstructionMask;

        // If using packages, transform mask to first stack (package) before cropping.
        if (_parameters.packages.size() > 0)
        {
            RealImage m = *reconstructionMask;
            RealImage t = _stacks[_parameters.templateNumber];
            if (_verbose)
                ReconstructionModel::SaveReconImageFile(t, "templateStack");
            TransformMask(t, m, RigidTransformation());
            RealImage mt = t * m;
            CropImage(mt, m);
            maskCropped = mt;
            if (_verbose)
                ReconstructionModel::SaveReconImageFile(m, "templateMask");
        }
        else if (_parameters.structural || _stacks[_parameters.templateNumber].GetT() == 1) {
            // Original 3D behaviour
            RealImage m = *reconstructionMask;
            RealImage t = _stacks[_parameters.templateNumber];
            if (_verbose) ReconstructionModel::SaveReconImageFile(t, "templateStack");
            if (_verbose) ReconstructionModel::SaveReconImageFile(m, "templateMask");
            TransformMask(t, m, RigidTransformation());
            RealImage mt = t * m;
            CropImage(mt, m);
            maskCropped = mt;
        }
        else {
            // Original 4D behaviour
            CropImage(maskCropped, *reconstructionMask);
        }
        
        if (_verbose) ReconstructionModel::SaveReconImageFile(*reconstructionMask,"mask"); 
        if (_verbose) ReconstructionModel::SaveReconImageFile(maskCropped,"maskCropped"); 
    }


    void ReconstructionAlgorithm::CalculateSliceThickness()
    {
        stringstream ss;
        ss << "Slice thickness is:";
        if (_parameters.thickness.empty())
        {
            for (size_t i = 0; i < _stacks.size(); i++)
            {
                double dx, dy, dz;
                _stacks[i].GetPixelSize(&dx, &dy, &dz);
                _parameters.thickness.push_back(dz);
                ss << " " << _parameters.thickness[i];
            }
        }
        else
        {
            for (size_t i = 0; i < _stacks.size(); i++)
                ss << _parameters.thickness[i] << " ";
        }
        logger::debug(ss.str());
    }
    void ReconstructionAlgorithm::SplitStacksIntoPackages()
    {
        logger::debug("Splitting stacks into packages ... ");

        Array<double> newThickness;
        Array<int> newPackages;
        Array<RigidTransformation> newStackTransformations;
        Array<RealImage> newStacks;

        int q = 0;
        for (size_t i = 0; i < _stacks.size(); i++)
        {
            Array<RealImage> out_packages;
            if (_parameters.packages[i] > 1)
            {
                SplitImage(_stacks[i], _parameters.packages[i], out_packages);
                for (size_t j = 0; j < out_packages.size(); j++)
                {
                    newStacks.push_back(out_packages[j]);
                    q = q + 1;
                    if (_parameters.thickness.size() > 0)
                        newThickness.push_back(_parameters.thickness[i]);
                    newPackages.push_back(1);
                    if (!_stackTransformations.empty())
                        newStackTransformations.push_back(_stackTransformations[i]);
                }
            }
            else
            {
                newStacks.push_back(_stacks[i]);
                if (_parameters.thickness.size() > 0)
                    newThickness.push_back(_parameters.thickness[i]);
                newPackages.push_back(1);
                if (!_stackTransformations.empty())
                    newStackTransformations.push_back(_stackTransformations[i]);
            }
        }

        _stacks = move(newStacks);
        _parameters.thickness = move(newThickness);
        _parameters.packages = move(newPackages);
        _stackTransformations = move(newStackTransformations);

        if (_verbose)
        {
            for (int i = 0; i < _stacks.size(); ++i)
            {
                _stacks[i].Write((boost::format("stack%1%_singlepackage.nii.gz") % i).str().c_str());
            }
        }

        logger::debug((boost::format("New number of stacks : %1%") % _stacks.size()).str());
    }

    void ReconstructionAlgorithm::CalculateDisplacements(int &iter, RealImage& reconstructed)
    {
        logger::trace("Reconstruction4DFunction::CalculateDisplacements");
        // Calculate Displacement
        _meanDisplacement.push_back(_reconstruction.CalculateDisplacement());
        _meanWeightedDisplacement.push_back(_reconstruction.CalculateWeightedDisplacement());
    }

    void ReconstructionAlgorithm::SaveIntermediateSimulatedSlices(int iter, int i)
    {
        // Save intermediate simulated slices
        if (_verbose)
        {
            if (_parameters.intensityMatchingMode >= IM_MODE::FULL_RESCALE_EXCEPT_STACKS)
            {
                _reconstruction.CalculateCorrectedSlices();
                _reconstruction.SaveCorrectedSlicesAsStack(_stacks, iter, i );
                if (_parameters.sigma > 0)
                    _reconstruction.SaveBiasFieldsAsStack(_stacks, iter, i );
            }
            _reconstruction.SaveSimulatedSlicesAsStack(_stacks, iter, i );
            _reconstruction.CalculateError();
            _reconstruction.SaveErrorAsStack(_stacks, iter, i);
        }
    }
    
    void ReconstructionAlgorithm::SaveIntermediateSimulatedSlices(int iter)
    {
        if (_verbose)
        {
            _reconstruction.SaveSimulatedSlicesAsStack(_stacks, iter);
            _reconstruction.SaveSimulatedWeights(_stacks, iter);
            _reconstruction.CalculateError();
            _reconstruction.SaveErrorAsStack(_stacks, iter);
        }
    }

    int ReconstructionAlgorithm::CalcReconstuctionIterations(int iter)
    {
        // if last iteration
        if (iter == _parameters.iterations - 1)
            return _parameters.srIterations;
        else if (_parameters.srIterations <= 0)
            return 0;
        else if (_parameters.srIterations <= 2)
            return 1;
        else
            return floor(_parameters.srIterations / 3);
    }

    //amount of smoothing (given by lambda) is decreased with improving alignment
    //delta (to determine edges) stays constant throughout
    void ReconstructionAlgorithm::SetSmoothingParameters(int iter)
    {
        if (iter == _parameters.iterations - 1)
            _reconstruction.SetSmoothingParameters(_parameters.delta, _parameters.lastIterLambda);
        else
        {
            double l = _parameters.lambda;
            for (int i = 0; i < _parameters.levels; i++)
            {
                if (iter == _parameters.iterations * (_parameters.levels - i - 1) / _parameters.levels)
                    _reconstruction.SetSmoothingParameters(_parameters.delta, l);
                l *= 2;
            }
        }
    }

    void ReconstructionAlgorithm::SliceToVolumeRegistration(int iter)
    {
        logger::trace("SliceToVolumeRegistration");
        if (iter > 0)
        {
            logger::trace((boost::format("Iteration %1%") % iter).str());

            // Determine if SVR masking option is enabled.
            bool svrMaskSlicesOption = _parameters.sliceMaskingMode == MASK_MODE::MASK_SLICES_ONLY_SVR || 
                                 _parameters.sliceMaskingMode == MASK_MODE::MASK_SLICES_ALWAYS;

            _reconstruction.SliceToVolumeRegistration(iter, svrMaskSlicesOption);
        }
    }

}
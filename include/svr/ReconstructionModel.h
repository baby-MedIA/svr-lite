// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/Common.h"

using namespace std;
using namespace mirtk;
using namespace svr::Utility;

namespace svr {

    // Forward declarations
    namespace Parallel {
        class StackRegistrations;
        class CoeffInit;
        class Superresolution;
        class MStep;
        class EStep;
        class SStep;
        class Bias;
        class Scale;
        class NormaliseBias;
        class SimulateSlices;
        class Average;
        class SliceToVolumeRegistration;
        class SimulateStacks;
        class CalculateError;
        class CalculateCorrectedSlices;
        class AdaptiveRegularizationStepA;
        class AdaptiveRegularizationStepB;
    }

    /**
     * @brief Reconstruction abstract base class
     */
    class ReconstructionModel {
    protected:

        /********
         * Core *
         ********/
        RECON_TYPE _recon_type; // Reconstruction type
        bool _debug = false;            // Debug mode
        bool _profile = false;          // Profile output mode
        bool _verbose = false;          // Verbose mode

        /**********
         * Slices *
         **********/
        Array<RealImage> _slices;
        Array<RealImage> _simulated_slices;
        Array<RealImage> _simulated_weights;
        Array<RealImage> _simulated_inside;
        Array<RealImage> _slice_masks;
        Array<RealImage> _not_masked_slices;
        Array<RealImage> _error;
        Array<RealImage> _corrected_slices;
        Array<double> _slice_dt;        // Slice Temporal Resolution

        Array<int> _small_slices;       // slices identify as too small to be used
        deque<bool> _slice_inside;      // Indicator whether slice has an overlap with volumetric mask
        deque<bool> _slice_excluded;    // Slice Excluded
        Array<int> _force_excluded;     // forced excluded slices
        Array<int> _force_excluded_stacks; // Force exclusion of stacks
        Array<double> _slice_tre;       // Target Registration Error
        Array<int> _zero_slices;

        Array<double> _stack_factor;    // to restore original signal intensity of the MRI slices
        Array<int> _stack_index;        // Index of the stack from which the slice originated
        Array<int> _stack_loc_index;    // index of slice location in stack
        Array<int> _loc_index;          // running index of all 2D slice locations

        // Displacement
        Array<double> _slice_displacement;
        Array<double> _slice_tx;
        Array<double> _slice_ty;
        Array<double> _slice_tz;
        Array<double> _slice_weighted_displacement;

        /*************************
         * Slice Transformations *
         *************************/
        Array<RigidTransformation> _transformations;
        Array<RigidTransformation> _transformationsRwithMB;
        Array<RigidTransformation> _ref_transformations;

        /*************************
         * Reconstructed Volumes *
         *************************/
        RealImage _reconstructed;
        RealImage _reconstructed_mask;
        RealImage _evaluation_mask; //Future improvement: retire in favour of _mask.
        RealImage _brain_probability;
        Array<RealImage> _probability_maps;

        /***********************
         * Volume Coefficients *
         ***********************/
        Array<SLICECOEFFS> _volcoeffs;  // Structures to store the matrix of transformation between volume and slices

        /***********
         * Weights *
         ***********/
        RealImage _volume_weights;      // Weights for Gaussian reconstruction
        RealImage _confidence_map;      // Weights for regularization
        Array<RealImage> _weights;      // Voxel posteriors
        Array<double> _slice_weight;    // Slice posteriors
        
        /*********
         * Flags *
         *********/
        bool _template_flag;
        bool _no_sr;
        bool _masked_stacks;
        bool _no_masking_background;
        bool _template_created;     // Flag to say whether the template volume has been created
        bool _have_mask;            // Flag to say whether we have a mask
        bool _adaptive;             // use adaptive or non-adaptive regularisation (default:false)
        bool _robust_slices_only;   //do not exclude voxels, only whole slices
        bool _no_ts;

        /*************
         * Variables *
         *************/
        int _nmi_bins;          /// NMI registration bins
        double _quality_factor; /// Quality factor - higher means slower and better
        double _average_value;

        /***************************
         * EM Algorithm Parameters *
         ***************************/
        double _sigma;      // Variance for inlier voxel errors    
        double _mix;        // Proportion of inlier voxels        
        double _m;          // Uniform distribution for outlier voxels        
        double _mean_s;     // Mean for inlier slice errors        
        double _sigma_s;    // Variance for inlier slice errors        
        double _mean_s2;    // Mean for outlier slice errors        
        double _sigma_s2;   // Variance for outlier slice errors        
        double _mix_s;      // Proportion of inlier slices        
        double _step;       // Step size for likelihood calculation   
        
        /// Intensity min and max
        double _max_intensity;
        double _min_intensity;

        /******************
         * Bias and Scale *
         ******************/
        double _sigma_bias;     // Variance for bias field
        Array<RealImage> _bias; // Slice-dependent bias fields
        Array<double> _scale;   // Slice-dependent scales
        bool _global_bias_correction;// global bias correction flag
        double _low_intensity_cutoff;// low intensity cutoff for bias field estimation

        /**************************************************
         * Gradient descent and regularisation parameters *
         **************************************************/
        double _alpha;  // Step for gradient descent
        double _delta;  // Determine what is en edge in edge-preserving smoothing
        double _lambda; // Amount of smoothing
        double _average_volume_weight;  // Average voxel wights to modulate parameter alpha


        /************************
         * Structural exclusion *
         ************************/
        bool _structural;                   // perform structural based exclusion of slices and voxels
        int _current_iteration;
        double _global_NCC_threshold;
        Array<int> _structural_slice_weight;
        Array<RealImage> _slice_ssim_maps;
        int _local_SSIM_window_size;
        double _local_SSIM_threshold;

        /********************
         * Voxel Neighbours *
         ********************/
        const int _directions[13][3] = {
            {1, 0, -1},
            {0, 1, -1},
            {1, 1, -1},
            {1, -1, -1},
            {1, 0, 0},
            {0, 1, 0},
            {1, 1, 0},
            {1, -1, 0},
            {1, 0, 1},
            {0, 1, 1},
            {1, 1, 1},
            {1, -1, 1},
            {0, 0, 1}
        };

        /*********************
         * Parallel Functors *
         *********************/
        friend class Parallel::StackRegistrations;
        friend class Parallel::CoeffInit;
        friend class Parallel::Superresolution;
        friend class Parallel::MStep;
        friend class Parallel::EStep;
        friend class Parallel::SStep;
        friend class Parallel::Bias;
        friend class Parallel::Scale;
        friend class Parallel::NormaliseBias;
        friend class Parallel::SimulateSlices;
        friend class Parallel::Average;
        friend class Parallel::SliceToVolumeRegistration;
        friend class Parallel::SimulateStacks;
        friend class Parallel::NormaliseBias;
        friend class Parallel::Superresolution;
        friend class Parallel::CalculateError;
        friend class Parallel::CalculateCorrectedSlices;
        friend class Parallel::AdaptiveRegularizationStepA;
        friend class Parallel::AdaptiveRegularizationStepB;

        /// Read Transformations common function
        void ReadTransformations(const char *folder, size_t file_count, Array<RigidTransformation>& transformations);

        // Common calculate function
        double Calculate(const RigidTransformation& drift, double (ReconstructionModel::*)());

        // Tukey Window Function
        inline double wintukey(double angdiff, double alpha) {
            if (fabs(angdiff) > PI * (1 - alpha))
                return (1 + cos((fabs(angdiff) - PI * (1 - alpha)) / alpha)) / 2;
            return 1;
        }

        // Register stacks in 3D-only
        void StackRegistrations3D(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations, RealImage& target);
        
        // Get 3D version of 4D stacks
        void Get3DStacks(const Array<RealImage>& stacks, Array<RealImage>& sim3dStacks);

    public:

        /// Reconstruction constructor
        ReconstructionModel();
        
        /// Reconstruction destructor
        virtual ~ReconstructionModel(); 

        int _number_of_slices_org;
        double _average_thickness_org;
        
        /**
         * @brief Remember volumetric mask and smooth it if necessary.
         * @param mask Mask to be set.
         * @param sigma
         * @param threshold
         * @return The resulting resolution of the template image.
         */
        void SetMask(RealImage& mask, double sigma, double threshold = 0.5);

        /// Mask slices based on the reconstruction mask
        void MaskSlices();

        // Create list of slice masks
        void CreateSliceMasks();
        
        /**
         * @brief Create average image from the stacks and volumetric transformations.
         * @param stacks
         * @param stack_transformations
         * @return
         */
        RealImage CreateAverage(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations);

        /**
         * @brief Calculate initial registrations.
         * @param stacks
         * @param stack_transformations
         * @param templateNumber
         */
        void StackRegistrations(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations);
        void StackRegistrations(Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations, const int templateNumber);
        void StackRegistrationsWithTemplate(const Array<RealImage>& stacks, Array<RigidTransformation>& stack_transformations, RealImage input_template);


        /**
         * @brief Saves probability map.
         * @param i
         */
        void SaveProbabilityMap(int i);

        /// Initialise values of variables and parameters for EM
        void InitializeEMValues();

        /// Perform EStep for calculation of voxel-wise and slice-wise posteriors (weights)
        void EStep();

        /// Calculate slice-dependent scale
        void Scale();

        // Calculate slice-dependent bias fields
        /// Run slice bias correction
        void Bias();

        /// Run MStep
        void MStep(int iter);
        

        /// ***********************
        /// Quality Evaluation Functions
        /// ***********************
        double CalculateSliceError(const int& index);
        double Consistency(bool exclude_slices); // SSD
        double CalculateRSE(bool exclude_slices);
        double CalculateNRMSE(bool exclude_slices);
        double CalculateNCC(bool exclude_slices);
        void SliceSummary(int& included_count, int& excluded_count, int& outside_count);

        /// ***********************

        /// Read transformations
        void ReadTransformations(const char *folder);

        /// Restore slice intensities to their original values to recover original scaling
        void RestoreSliceIntensities();

        // Compute local SSIM map between simulated and original slices
        void SliceSSIMMap(int inputIndex);

        void ReadSliceTransformations(const char *folder);

        /// Determine resolution of volume to reconstruct
        double GetReconstructedResolutionFromTemplateStack(const RealImage& stack);

        /**
         * @brief Initialise reconstructed volume from static mask.
         * @param mask
         * @param resolution
         * @param temporalSize size dimension
         * @param temporalResolution
         */
        void CreateTemplateFromStaticMask(const RealImage& mask, double resolution, int temporalSize, double temporalResolution);

        
        void CreateTemplate(ImageAttributes attr, 
                            const RealImage& templateImage,
                            double resolution,
                            int temporalSize,
                            double temporalResolution);


        /**
         * @brief Match stack intensities with respect to the masked ROI.
         * @param stacks
         * @param stack_transformations
         * @param averageValue
         * @param together
         */
        void MatchStackIntensitiesWithMasking(Array<RealImage>& stacks,
            const Array<RigidTransformation>& stack_transformations,
            double averageValue,
            bool together = false);

        /**
         * @brief Create slices from the stacks and slice-dependent transformations from stack transformations
         * @param stacks
         * @param stack_transformations
         * @param thickness
         * @param probability_maps
         */
        void CreateSlicesAndTransformations(const Array<RealImage>& stacks,
            const Array<RigidTransformation>& stack_transformations,
            const Array<double>& thickness,
            const Array<RealImage>& probability_maps = Array<RealImage>());

        /// Calculate transformation matrix between slices and voxels
        void CoeffInit();

        /// PSF-weighted reconstruction
        void GaussianReconstruction();

        /// Initialise variables and parameters for EM
        virtual void InitializeEM();
        
        /// Initalise parameters of EM robust statistics
        virtual void InitializeRobustStatistics();

        /// Calculate corrected slices
        void CalculateCorrectedSlices();

        /// Simulate slices
        void SimulateSlices();

        /// Simulate stacks
        void SimulateStacks();

        /// Normalise bias
        void NormaliseBias(int iter, int rec_inter);

        /// Slice-to-volume registration
        void SliceToVolumeRegistration(int iter = 0, bool maskSlices = false);
        
        /**
         * @brief Volume-to-volume Registration
         * @param target
         * @param source
         * @param rigidTransf
         */
        void VolumeToVolumeRegistration(const GreyImage& target, const GreyImage& source, RigidTransformation& rigidTransf);

        /// Calculate displacement
        double CalculateDisplacement();
        /// Calculate displacement relative to alignment
        inline double CalculateDisplacement(const RigidTransformation& drift) { return Calculate(drift, &ReconstructionModel::CalculateDisplacement); }
        /// Calculate weighted displacement
        double CalculateWeightedDisplacement();
        /// Calculate weighted displacement relative to alignment
        inline double CalculateWeightedDisplacement(const RigidTransformation& drift) { return Calculate(drift, &ReconstructionModel::CalculateWeightedDisplacement); }

        /// Initialise target registration error (TRE)
        inline void InitTRE() { _slice_tre = Array<double>(_slices.size(), -1); }
        /// Calculate target registration error (TRE)
        double CalculateTRE();
        /// Calculate target registration error (TRE) relative to alignment
        inline double CalculateTRE(const RigidTransformation& drift) { return Calculate(drift, &ReconstructionModel::CalculateTRE); }

        /**
         * @brief Smooth transformations.
         * @param sigma
         * @param niter
         * @param use_slice_inside
         */
        void SmoothTransformations(double sigma, int niter = 10, bool use_slice_inside = false);

        /// Scale transformations
        void ScaleTransformations(double scale);

        /// Apply static mask to reconstructed 4d volume
        inline void StaticMaskReconstructedVolume() { StaticMaskVolume(_reconstructed, _reconstructed_mask, -1); }

        /// Calculate error
        void CalculateError();

        /// Run superresolution reconstruction step
        void Superresolution(int iter);

        /// Edge-preserving regularization with confidence map
        void AdaptiveRegularization(int iter, const RealImage& original);

        /// Read reference transformations
        void ReadRefTransformations(const char *folder);

        /// Calculate entropy
       double CalculateEntropy();

        /// Run structure-based rejection of outliers
        void GlobalStructuralExclusion();

        double LocalSSIM(const RealImage& slice, const RealImage& sim_slice, const RealImage& mask, const int t);

        /// Run SStep 
        void SStep();


        ////////////////////////////////////////////////////////////////////////////////
        // Inline/template definitions
        ////////////////////////////////////////////////////////////////////////////////

        /// Zero-mean Gaussian PDF
        inline double G(double x, double s) {
            return _step * exp(-x * x / (2 * s)) / (sqrt(6.28 * s));
        }

        /// Uniform PDF
        inline double M(double m) {
            return m * _step;
        }
                
        /// Return mask
        inline const RealImage& GetMask() {
            return _reconstructed_mask;
        }

        /// Return volume weights
        inline const RealImage& GetVolumeWeights() {
            return _volume_weights;
        }

        /// Enable debug mode
        inline void DebugOn() {
            _debug = true;
        }

        /// Disable debug mode
        inline void DebugOff() {
            _debug = false;
            _verbose = false;
        }

        /// Enable profiling mode
        inline void ProfileOn() {
            logger::init_profile_log();
            _profile = true;
        }

        /// Disable profiling mode
        inline void ProfileOff() {
            _profile = false;
        }

        /// Enable verbose mode
        //(all logging controlled by logger class, but this can be used for efficiency gains) 
        inline void VerboseOn() {
            _verbose = true;
            ProfileOn();
        }

        /// Set blurring flag
        // inline void SetBlurring(bool flag_blurring) {
        //     _blurring = flag_blurring;
        // }
        
               
        /// Set sigma flag
        inline void SetSigma(double sigma) {
            _sigma_bias = sigma;
        }

        /// Use faster lower quality reconstruction
        inline void SpeedupOn() {
            _quality_factor = 1;
        }

        /// Use slower better quality reconstruction
        inline void SpeedupOff() {
            _quality_factor = 2;
        }

        /// Enable global bias correction
        inline void GlobalBiasCorrectionOn() {
            _global_bias_correction = true;
        }

        /// Disable global bias correction
        inline void GlobalBiasCorrectionOff() {
            _global_bias_correction = false;
        }

        /// Enable robust slices mode
        inline void ExcludeWholeSlicesOnly() {
            _robust_slices_only = true;
            cout << "Exclude only whole slices." << endl;
        }


        inline void SetCurrentIteration(int current_iteration) {
            _current_iteration = current_iteration;
        }

        /// Set structural flag
        inline void SetStructural(bool flag_structural) {
            _structural = flag_structural;
        }

        /// Set global NCC used for structural exclusion
        inline void SetGlobalNCC(double global_NCC_threshold) {
            _global_NCC_threshold = global_NCC_threshold;
        }

        /// Set smoothing parameters
        inline void SetSmoothingParameters(double delta, double lambda) {
            _delta = delta;
            _lambda = lambda * delta * delta;
            _alpha = 0.05 / lambda;
            if (_alpha > 1)
                _alpha = 1;

            logger::trace((boost::format("delta = %1%  lambda = %2% alpha = %3%") % _delta % lambda % _alpha).str());
        }

        /// Set slices which need to be excluded by default
        inline void SetForceExcludedSlices(const Array<int>& force_excluded) {
            _force_excluded = force_excluded;
        }

        /// Set force excluded stacks
        inline void SetForceExcludedStacks(const Array<int>& force_excluded_stacks) {
            _force_excluded_stacks = force_excluded_stacks;
        }

        inline bool InSliceMask(const int& sliceIdx, const int& x, const int& y) {
            return _slice_masks[sliceIdx](x,y,0) == 1;
        }

        /// Return reconstructed 4d volume
        inline const RealImage& GetReconstructed() {
            return _reconstructed;
        }
        
        /// Set reconstructed volume
        inline void SetReconstructed(const RealImage& reconstructed) {
            _reconstructed = reconstructed;
            _template_created = true;
        }

        /// Initialise stack factor
        inline void InitStackFactor(const Array<RealImage>& stacks) {
            _stack_factor = Array<double>(stacks.size(), 1);
        }

        inline void OutputStackFactor(){
            cout << "Stack factors: ";
            for (auto sf : _stack_factor){
                cout << sf << ", ";
            }
            cout << endl;
        }

        /// Reset slices and transformations
        void ResetSlicesAndTransformations() {
            _slices.clear();
            _simulated_slices.clear();
            _simulated_weights.clear();
            _simulated_inside.clear();
            _stack_index.clear();
            _loc_index.clear();
            _stack_loc_index.clear();
            _transformations.clear();
            _slice_excluded.clear();
            _probability_maps.clear();
        }

        /// Initialise corrected slices
        inline void InitCorrectedSlices() {
            _corrected_slices = _slices;
        }

        /// Initialise error
        inline void InitError() {
            _error = _slices;
        }


        /// Scale volume to match the slice intensities
        inline void ScaleVolume() {
            logger::trace("Scaling volume: ");

            ScaleVolume(_reconstructed);
        }

        /// Scale volume common function
        virtual void ScaleVolume(RealImage& reconstructed);

        inline vector<double> GetSliceScale(){
            return _scale;
        }

        /// ***********************
        /// Save functions
        /// ***********************

        /// Core save functions that use common naming convention
        static void SaveStackFile(const RealImage& image, 
                              const string& suffix, 
                              const int stackNum = -1,
                              const int iterNum = -1,
                              const int recIterNum = -1);
        static void SaveReconImageFile( const RealImage& image, 
                                    const string& suffix, 
                                    const int iterNum = -1,
                                    const int recIterNum = -1);
        static void SaveReconImageFile( const RealImage& image, 
                                    const int iterNum,
                                    const int recIterNum = -1);
        static void SaveSliceFile(const RealImage& image, const string& suffix, const int sliceNum);
        void SaveStacks(const Array<RealImage>& save_stacks, const string& name);
        void SaveAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, const string& name, int iter = -1, int rec_iter = -1);
        void Save2dAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, const string& name, int iter = -1, int rec_iter = -1);
        void GetAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, Array<RealImage>& output);
        void Get2dAsStack(const Array<RealImage>& source, const Array<RealImage>& stacks, Array<RealImage>& output);

        /// Save slices
        void SaveSlices();
        void SaveSlicesAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);
        
        /// Save simulated slices
        void SaveSimulatedSlices();
        void SaveSimulatedSlicesAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);
        void SaveScaledSimulatedSlicesAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);
        void GetSimulatedSlicesAsStacks(const vector<RealImage>& stacks, vector<RealImage>& simulatedStacks);
        void GetScaledSimulatedSlicesAsStacks(const vector<RealImage>& stacks, vector<RealImage>& simulatedStacks);

        /// Save simulated slices
        void SaveSliceMasks();

        /// Save weights
        void SaveWeights();
        void SaveWeightsAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);

        /// Save transformations
        void SaveTransformations();
        static void SaveStackTransformation(const RigidTransformation& transform,const int stackNum);

        /// Save simulated weights
        void SaveSimulatedWeights(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);


        /// Save corrected slices
        void SaveCorrectedSlices();
        void SaveCorrectedSlicesAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);

        /// Save error
        void SaveError();
        void SaveErrorAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1);

        /// Save bias fields
        void SaveBiasFields();
        void SaveBiasFieldsAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1) ;
        void GetBiasFieldsAsStack(const vector<RealImage>& stacks, vector<RealImage>& output);

        /// Save ssim maps
        void SaveSSIMMapsAsStack(const Array<RealImage>& stacks, int iter = -1, int rec_iter = -1) ;

        /// Save reconstruction metadata as JSON file
        void SaveMetaData(const boost::timer::cpu_timer& reconstruction_timer);

        /// Save slice info
        void SaveSliceInfo(const char *filename, const Array<string>& stack_filenames);

    };  // end of Reconstruction class definition
    
} // namespace svr

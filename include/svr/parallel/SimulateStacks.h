// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/ReconstructionModel.h"

using namespace mirtk;
using namespace svr;
using namespace svr::Utility;

namespace svr::Parallel {

    // //-------------------------------------------------------------------
    class SimulateStacks {
        ReconstructionModel *reconstructor;

    public:
        SimulateStacks(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            RealImage PSF, tPSF;
            //get resolution of the volume
            double vx, vy, vz;
            reconstructor->_reconstructed.GetPixelSize(&vx, &vy, &vz);
            //volume is always isotropic
            const double& res = vx;

            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                if (reconstructor->_slice_excluded[inputIndex])
                    continue;

                //read the slice
                RealImage& slice = reconstructor->_slices[inputIndex];

                //prepare structures for storage
                SLICECOEFFS slicecoeffs(slice.GetX(), Array<VOXELCOEFFS>(slice.GetY()));

                //PSF will be calculated in slice space in higher resolution

                //get slice voxel size to define PSF
                double dx, dy, dz;
                slice.GetPixelSize(&dx, &dy, &dz);

                //sigma of 3D Gaussian (sinc with FWHM=dx or dy in-plane, Gaussian with FWHM = dz through-plane)

                double sigmax, sigmay, sigmaz;
                if (reconstructor->_recon_type == _3D) {
                    sigmax = 1.2 * dx / 2.3548;
                    sigmay = 1.2 * dy / 2.3548;
                    sigmaz = dz / 2.3548;
                }

                if (reconstructor->_recon_type == _1D) {
                    sigmax = 0.5 * dx / 2.3548;
                    sigmay = 0.5 * dy / 2.3548;
                    sigmaz = dz / 2.3548;
                }

                if (reconstructor->_recon_type == _interpolate) {
                    sigmax = 0.5 * dx / 2.3548;
                    sigmay = 0.5 * dx / 2.3548;
                    sigmaz = 0.5 * dx / 2.3548;
                }

                //calculate discretized PSF

                //isotropic voxel size of PSF - derived from resolution of reconstructed volume
                double size = res / reconstructor->_quality_factor;

                //number of voxels in each direction
                //the ROI is 2*voxel dimension

                int xDim = round(2 * dx / size);
                int yDim = round(2 * dy / size);
                int zDim = round(2 * dz / size);
                ///test to make dimension always odd
                xDim = xDim / 2 * 2 + 1;
                yDim = yDim / 2 * 2 + 1;
                zDim = zDim / 2 * 2 + 1;
                ///end test

                //image corresponding to PSF
                ImageAttributes attr;
                attr._x = xDim;
                attr._y = yDim;
                attr._z = zDim;
                attr._dx = size;
                attr._dy = size;
                attr._dz = size;
                PSF.Initialize(attr);

                //centre of PSF
                double cx, cy, cz;
                cx = 0.5 * (xDim - 1);
                cy = 0.5 * (yDim - 1);
                cz = 0.5 * (zDim - 1);
                PSF.ImageToWorld(cx, cy, cz);

                double x, y, z;
                double sum = 0;
                int i, j, k;
                for (i = 0; i < xDim; i++)
                    for (j = 0; j < yDim; j++)
                        for (k = 0; k < zDim; k++) {
                            x = i;
                            y = j;
                            z = k;
                            PSF.ImageToWorld(x, y, z);
                            x -= cx;
                            y -= cy;
                            z -= cz;
                            //continuous PSF does not need to be normalized as discrete will be
                            PSF(i, j, k) = exp(-x * x / (2 * sigmax * sigmax) - y * y / (2 * sigmay * sigmay) - z * z / (2 * sigmaz * sigmaz));
                            sum += PSF(i, j, k);
                        }
                PSF /= sum;

                if (reconstructor->_verbose)
                    if (inputIndex == 0)
                        PSF.Write("slice0-PSF.nii.gz");

                //prepare storage for PSF transformed and resampled to the space of reconstructed volume
                //maximum dim of rotated kernel - the next higher odd integer plus two to accound for rounding error of tx,ty,tz.
                //Note conversion from PSF image coordinates to tPSF image coordinates *size/res
                int dim = (floor(ceil(sqrt(double(xDim * xDim + yDim * yDim + zDim * zDim)) * size / res) / 2)) * 2 + 1 + 2;
                //prepare image attributes. Voxel dimension will be taken from the reconstructed volume
                attr._x = dim;
                attr._y = dim;
                attr._z = dim;
                attr._dx = res;
                attr._dy = res;
                attr._dz = res;
                //create matrix from transformed PSF
                tPSF.Initialize(attr);
                //calculate centre of tPSF in image coordinates
                int centre = (dim - 1) / 2;

                //for each voxel in current slice calculate matrix coefficients
                int ii, jj, kk;
                int tx, ty, tz;
                int nx, ny, nz;
                int l, m, n;
                for (i = 0; i < slice.GetX(); i++)
                    for (j = 0; j < slice.GetY(); j++)
                        if (reconstructor->InSliceMask(inputIndex,i, j)) {
                            //calculate centrepoint of slice voxel in volume space (tx,ty,tz)
                            x = i;
                            y = j;
                            z = 0;
                            slice.ImageToWorld(x, y, z);
                            reconstructor->_transformations[inputIndex].Transform(x, y, z);
                            reconstructor->_reconstructed.WorldToImage(x, y, z);
                            tx = round(x);
                            ty = round(y);
                            tz = round(z);

                            //Clear the transformed PSF
                            memset(tPSF.Data(), 0, sizeof(RealPixel) * tPSF.NumberOfVoxels());

                            //for each POINT3D of the PSF
                            for (ii = 0; ii < xDim; ii++)
                                for (jj = 0; jj < yDim; jj++)
                                    for (kk = 0; kk < zDim; kk++) {
                                        //Calculate the position of the POINT3D of
                                        //PSF centred over current slice voxel
                                        //This is a bit complicated because slices
                                        //can be oriented in any direction

                                        //PSF image coordinates
                                        x = ii;
                                        y = jj;
                                        z = kk;
                                        //change to PSF world coordinates - now real sizes in mm
                                        PSF.ImageToWorld(x, y, z);
                                        //centre around the centrepoint of the PSF
                                        x -= cx;
                                        y -= cy;
                                        z -= cz;

                                        //Need to convert (x,y,z) to slice image
                                        //coordinates because slices can have
                                        //transformations included in them (they are
                                        //nifti)  and those are not reflected in
                                        //PSF. In slice image coordinates we are
                                        //sure that z is through-plane

                                        //adjust according to voxel size
                                        x /= dx;
                                        y /= dy;
                                        z /= dz;
                                        //center over current voxel
                                        x += i;
                                        y += j;

                                        //convert from slice image coordinates to world coordinates
                                        slice.ImageToWorld(x, y, z);

                                        //x+=(vx-cx); y+=(vy-cy); z+=(vz-cz);
                                        //Transform to space of reconstructed volume
                                        reconstructor->_transformations[inputIndex].Transform(x, y, z);
                                        //Change to image coordinates
                                        reconstructor->_reconstructed.WorldToImage(x, y, z);

                                        //determine coefficients of volume voxels for position x,y,z
                                        //using linear interpolation

                                        //Find the 8 closest volume voxels

                                        //lowest corner of the cube
                                        nx = (int)floor(x);
                                        ny = (int)floor(y);
                                        nz = (int)floor(z);

                                        //not all neighbours might be in ROI, thus we need to normalize
                                        //(l,m,n) are image coordinates of 8 neighbours in volume space
                                        //for each we check whether it is in volume
                                        sum = 0;
                                        //to find wether the current slice voxel has overlap with ROI
                                        bool inside = false;
                                        for (l = nx; l <= nx + 1; l++)
                                            if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
                                                for (m = ny; m <= ny + 1; m++)
                                                    if ((m >= 0) && (m < reconstructor->_reconstructed.GetY()))
                                                        for (n = nz; n <= nz + 1; n++)
                                                            if ((n >= 0) && (n < reconstructor->_reconstructed.GetZ())) {
                                                                const double weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z));
                                                                sum += weight;
                                                                if (reconstructor->_reconstructed_mask(l, m, n) == 1)
                                                                    inside = true;
                                                            }
                                        //if there were no voxels do nothing
                                        if (sum <= 0 || !inside)
                                            continue;
                                        //now calculate the transformed PSF
                                        for (l = nx; l <= nx + 1; l++)
                                            if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
                                                for (m = ny; m <= ny + 1; m++)
                                                    if ((m >= 0) && (m < reconstructor->_reconstructed.GetY()))
                                                        for (n = nz; n <= nz + 1; n++)
                                                            if ((n >= 0) && (n < reconstructor->_reconstructed.GetZ())) {
                                                                const double weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z));

                                                                //image coordinates in tPSF
                                                                //(centre,centre,centre) in tPSF is aligned with (tx,ty,tz)
                                                                int aa = l - tx + centre;
                                                                int bb = m - ty + centre;
                                                                int cc = n - tz + centre;

                                                                //resulting value
                                                                double value = PSF(ii, jj, kk) * weight / sum;

                                                                //Check that we are in tPSF
                                                                if ((aa < 0) || (aa >= dim) || (bb < 0) || (bb >= dim) || (cc < 0) || (cc >= dim)) {
                                                                    stringstream err;
                                                                    err << "Error while trying to populate tPSF. " << aa << " " << bb << " " << cc << endl;
                                                                    err << l << " " << m << " " << n << endl;
                                                                    err << tx << " " << ty << " " << tz << endl;
                                                                    err << centre << endl;
                                                                    if (reconstructor->_verbose) tPSF.Write("tPSF.nii.gz");
                                                                    throw runtime_error(err.str());
                                                                } else
                                                                    //update transformed PSF
                                                                    tPSF(aa, bb, cc) += value;
                                                            }

                                    } //end of the loop for PSF points

                            //store tPSF values
                            for (ii = 0; ii < dim; ii++)
                                for (jj = 0; jj < dim; jj++)
                                    for (kk = 0; kk < dim; kk++)
                                        if (tPSF(ii, jj, kk) > 0) {
                                            POINT3D p;
                                            p.x = ii + tx - centre;
                                            p.y = jj + ty - centre;
                                            p.z = kk + tz - centre;
                                            p.value = tPSF(ii, jj, kk);
                                            slicecoeffs[i][j].push_back(p);
                                        }

                        } //end of loop for slice voxels

                //Calculate simulated slice
                reconstructor->_simulated_slices[inputIndex].Initialize(reconstructor->_slices[inputIndex].Attributes());
                reconstructor->_simulated_weights[inputIndex].Initialize(reconstructor->_slices[inputIndex].Attributes());

                // candidate for voxelfunction
                for (int l = 0; l < reconstructor->_slices[inputIndex].GetT(); l++)
                    for (int j = 0; j < reconstructor->_slices[inputIndex].GetY(); j++)
                        for (int i = 0; i < reconstructor->_slices[inputIndex].GetX(); i++)
                            if (reconstructor->InSliceMask(inputIndex,i,j)) {
                                double weight = 0;
                                const size_t n = slicecoeffs[i][j].size();
                                for (size_t k = 0; k < n; k++) {
                                    const POINT3D& p = slicecoeffs[i][j][k];
                                    reconstructor->_simulated_slices[inputIndex](i, j, 0, l) += p.value * reconstructor->_reconstructed(p.x, p.y, p.z, l);
                                    weight += p.value;
                                }
                                if (weight > 0) {
                                    reconstructor->_simulated_slices[inputIndex](i, j, 0, l) /= weight;
                                    reconstructor->_simulated_weights[inputIndex](i, j, 0, l) = weight;
                                }
                            }
            }  //end of loop through the slices
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}
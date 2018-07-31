/*
 * Created by Mykyta Haranko
 *  (2018 DESY)
 *
 *  email:mykyta.haranko@desy.de
 */

#ifndef CBCHITRECOVERY_H
#define CBCHITRECOVERY_H 1

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// system includes <>
#include <string>
#include <map>

// ROOT includes <>
#include "TObject.h"
#include <TROOT.h>
#include <TMath.h>
#include "TVector3.h"
#include "TVector2.h"


namespace eutelescope
{

    //! CMS CBC Hit Recovery processor for Marlin (recovers real track hits from the given position on the dut)

    class CBCHitRecovery : public marlin::Processor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new CBCHitRecovery;
	    }

	    //! Default constructor
	    CBCHitRecovery ();

	    virtual void init ();

	    virtual void processRunHeader (LCRunHeader * run);

	    virtual void processEvent (LCEvent * evt);

	    virtual void check (LCEvent * evt);

	    void bookHistos();

	    void fillHistos();

	    void vector_set_length(double *& vec, double length);
	    double vector_get_length(const double *vec);
	    double cos_alpha(const double *vec1,const double *vec2);
	    double* hit_pos(const double track[],const double normal[],const double virtual_hit[], double distance);

	    virtual void end();

	    std::string _InputTrackCollectionName;

	    std::string _InputFitHitsCollectionName;

	    std::string _cbcInputCollectionName;

	    std::string _cbcFitHitOutputCollectionName;

	    std::string _cbcTracksOutputCollectionName;

	    int _cbcVirtualDUTId;

	    std::vector<int> _cbcRealDUTsVec;

	    double *_VirtualDutNormal_zplus;
	    double *_VirtualDutNormal_zminus;
	    double *_VirtualDutPos;

	    gear::SiPlanesParameters * _siPlanesParameters;
	    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

	    // map of the vectors
	    std::map < int, double* > _dutPosMap;
	    std::map < int, double* > _dutNormalMap;

            //! Alignment collection name
            std::vector<std::string> _alignmentCollectionName;

            //! Prealignment collection name
            std::vector<std::string> _pre_alignmentCollectionName;

            //! Original prealignment applied to the original hit collection
            std::string _originalPreAlignment;


            //! Telescope plane count
            int _nTelPlanes;

            bool _useOriginalPreAlignment;

            //! Flag if alignment is loaded
            bool _alignmentloaded;

            //! Flag if prealignment is loaded
            bool _pre_alignmentloaded;

            //! Flag if reference is loaded
            bool _referenceloaded;

            //! The variables to store the alignment from the files
            double * _dut_align_x_1;
            double * _dut_align_y_1;
            double * _dut_align_z_1;
            double * _dut_align_a_1;
            double * _dut_align_b_1;
            double * _dut_align_c_1;
            double * _dut_align_x_error_1;
            double * _dut_align_y_error_1;
            double * _dut_align_z_error_1;
            double * _dut_align_a_error_1;
            double * _dut_align_b_error_1;
            double * _dut_align_c_error_1;
            double * _dut_pre_align_x_1;
            double * _dut_pre_align_y_1;
            double * _dut_pre_align_z_1;
            double * _dut_pre_align_a_1;
            double * _dut_pre_align_b_1;
            double * _dut_pre_align_c_1;
            double _dut_original_pre_align_x_1;
            double _dut_original_pre_align_y_1;
            double _dut_original_pre_align_z_1;
            double _dut_original_pre_align_a_1;
            double _dut_original_pre_align_b_1;
            double _dut_original_pre_align_c_1;
            //! DUT positions from the reference hit collection
            double _x_refhit_1;
            double _y_refhit_1;
            double _z_refhit_1;
            double _a_refhit_1;
            double _b_refhit_1;
            double _c_refhit_1;

            //! The variables to store the alignment from the files
            double * _dut_align_x_2;
            double * _dut_align_y_2;
            double * _dut_align_z_2;
            double * _dut_align_a_2;
            double * _dut_align_b_2;
            double * _dut_align_c_2;
            double * _dut_align_x_error_2;
            double * _dut_align_y_error_2;
            double * _dut_align_z_error_2;
            double * _dut_align_a_error_2;
            double * _dut_align_b_error_2;
            double * _dut_align_c_error_2;
            double * _dut_pre_align_x_2;
            double * _dut_pre_align_y_2;
            double * _dut_pre_align_z_2;
            double * _dut_pre_align_a_2;
            double * _dut_pre_align_b_2;
            double * _dut_pre_align_c_2;
            double _dut_original_pre_align_x_2;
            double _dut_original_pre_align_y_2;
            double _dut_original_pre_align_z_2;
            double _dut_original_pre_align_a_2;
            double _dut_original_pre_align_b_2;
            double _dut_original_pre_align_c_2;
            //! DUT positions from the reference hit collection
            double _x_refhit_2;
            double _y_refhit_2;
            double _z_refhit_2;
            double _a_refhit_2;
            double _b_refhit_2;
            double _c_refhit_2;


            //! The pitch and pixel count of the DUT
            double _pitchx;
            double _pitchy;
            int _pixelx;
            int _pixely;

            //! Additional alignment and rotations from the GEAR file
            std::vector<float > _DUTalign1;
            std::vector<float > _DUTalign2;

            // de-align the hit
            TVector2 DoDeAlignment1(TVector3 hit_input);
            TVector2 DoDeAlignment2(TVector3 hit_input);


            //! The function to get the alignment from file
            void getAlignment (LCEvent * event);

            //! The function to get the prealignment from file
            void getPreAlignment (LCEvent * event);

            //! The reference hit collection
            std::string _referencecollectionname;

            //! The function to get the reference hit collection
            void getReference (LCEvent * event);

            //! round value
            int roundValue(double value) {
                int integer = int(floor(value));
                double res = value - integer;
                if (res >= 0.5) return integer+1;
                else return integer;
            }

	protected:

	    std::map < std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    };

    //! A global instance of the processor
    CBCHitRecovery gCBCHitRecovery;

}

#endif

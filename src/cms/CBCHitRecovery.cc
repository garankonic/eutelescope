/*
 * Created by Mykyta Haranko
 *  (2018 DESY)
 *
 *  email:mykyta.haranko@desy.de
 */

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <UTIL/LCTOOLS.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <glob.h>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>

// eutelescope includes ""
#include "EUTELESCOPE.h"
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelReferenceHit.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"



// processor include
#include "CBCHitRecovery.h"

#define PI 3.14159265

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


CBCHitRecovery::CBCHitRecovery ( ) : Processor ( "CBCHitRecovery" ),
_aidaHistoMap ( )
{

    _description = "CBCHitRecovery recovers the real hit positions from the virtual dut and track.";

    registerProcessorParameter ( "InputTrackCollectionName", "The name of the Tracks collection we want to read", _InputTrackCollectionName, string ( "tracks_input" ) );

    registerProcessorParameter ( "InputFitHitsCollectionName", "The name of the fithits collection we want to read", _InputFitHitsCollectionName, string ( "fithits_input" ) );
    
    registerProcessorParameter ( "CBCInputCollectionName", "The name of the CBC collection we want to read", _cbcInputCollectionName, string ( "cbc_input" ) );

    registerProcessorParameter ( "CBCFitHitOutputCollectionName", "The name of the CBC data collection we want to write", _cbcFitHitOutputCollectionName, string ( "cbc_fithits_out" ) );

    registerProcessorParameter ( "CBCTracksOutputCollectionName", "The name of the CBC data collection we want to write", _cbcTracksOutputCollectionName, string ( "cbc_tracks_out" ) );

    registerProcessorParameter ( "CBCVirtualDUTId", "The id of the virtual CBC DUT", _cbcVirtualDUTId, static_cast <int> (0) );
    
    registerProcessorParameter ( "CBCRealDUTsVec", "The ids of the real DUTs (has to be two)", _cbcRealDUTsVec, std::vector < int > () );

    registerOptionalParameter ("OriginalPreAlignment", "Any pre-Prealignment applied to the original hit collection.", _originalPreAlignment, std::string("prealignment"));

    registerOptionalParameter ("PreAlignmentCollectionName", "To reconstruct the local DUT coordinates, all prealignment collection names which have been applied must be entered here. The order should correspond to their application, i.e.: prealignment1 prealignment2 etc.", _pre_alignmentCollectionName, std::vector < std::string > () );

    registerOptionalParameter ("AlignmentCollectionName", "To reconstruct the local DUT coordinates, all alignment collection names which have been applied must be entered here. The order should correspond to their application, i.e.: alignment1 alignment2 etc.", _alignmentCollectionName, std::vector < std::string > ());


    registerOptionalParameter ("ReferenceCollectionName", "The name of the reference collection. This is needed to get the correct Z DUT position after alignments.", _referencecollectionname, std::string("referenceHit"));

    registerOptionalParameter ("UseOriginalPreAlignment", "Subtract any pre-Prealignment applied to the original hit collection?", _useOriginalPreAlignment, static_cast < bool > (true) );

}


void CBCHitRecovery::init ( )
{
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	printParameters ( );

	_siPlanesParameters  = const_cast < gear::SiPlanesParameters* > ( & ( Global::GEAR -> getSiPlanesParameters ( ) ) );
	_siPlanesLayerLayout = const_cast < gear::SiPlanesLayerLayout* > ( & ( _siPlanesParameters -> getSiPlanesLayerLayout ( ) ) );

        // Take all layers defined in GEAR geometry
        _nTelPlanes = _siPlanesLayerLayout->getNLayers();

	for ( int i = 0; i < _siPlanesParameters -> getSiPlanesNumber ( ); i++ )
	{
		int cPlaneID = _siPlanesLayerLayout -> getID (i);
		if ( std::find(_cbcRealDUTsVec.begin(), _cbcRealDUTsVec.end(), cPlaneID) != _cbcRealDUTsVec.end())
		{
			// get the angles
			double alpha = _siPlanesLayerLayout -> getLayerRotationZY ( i );
			double beta  = _siPlanesLayerLayout -> getLayerRotationZX ( i );

			// we need to calculate the normal vectors 
			double *normalvector = new double[3];
			normalvector[0] = sin ( beta * PI / 180.0 );
    			normalvector[1] = -sin ( alpha * PI / 180.0 );
    			normalvector[2] = cos ( alpha * PI / 180.0 ) * cos ( beta * PI / 180.0 );

			// put it to the map now
			_dutNormalMap.insert ( make_pair ( cPlaneID, normalvector ) );

			// we need to calculate the normal vectors 
			double *posvector = new double[3];
			posvector[0] = _siPlanesLayerLayout -> getLayerPositionX( i );
    			posvector[1] = _siPlanesLayerLayout -> getLayerPositionY( i );
    			posvector[2] = _siPlanesLayerLayout -> getLayerPositionZ( i ) + 0.5 *_siPlanesLayerLayout -> getSensitiveThickness( i );

			// put it to the map now
			_dutPosMap.insert ( make_pair ( cPlaneID, posvector ) );

			// test if it was saved properly
			double *pos_test = _dutPosMap.at(cPlaneID);
			streamlog_out ( MESSAGE2 ) << "The position vector of the DUT" << cPlaneID << " is: X = " << pos_test[0] << ", Y = " << pos_test[1] << ", Z = " << pos_test[2] << endl;
			double *normal_test = _dutNormalMap.at(cPlaneID);
			streamlog_out ( MESSAGE2 ) << "The normal vector of the DUT" << cPlaneID << " is: X = " << normal_test[0] << ", Y = " << normal_test[1] << ", Z = " << normal_test[2] << endl;
                }

                // fill the picthes
                if(_cbcVirtualDUTId == cPlaneID) {
                    // get the DUT parameters from gear
                    _pitchx = _siPlanesLayerLayout->getSensitivePitchX(i);
                    _pitchy = _siPlanesLayerLayout->getSensitivePitchY(i);
                    _pixelx = _siPlanesLayerLayout->getSensitiveNpixelX(i);
                    _pixely = _siPlanesLayerLayout->getSensitiveNpixelY(i);
                }
                // fill geometry
                if(cPlaneID == _cbcRealDUTsVec.at(0)) {
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerPositionX(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT X shift according to gear file: " << _DUTalign1.at(0) << endl;
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerPositionY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT Y shift according to gear file: " << _DUTalign1.at(1) << endl;
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerPositionZ(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT Z shift according to gear file: " << _DUTalign1.at(2) << endl;
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerRotationZY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT ZY rotation according to gear file: " << _DUTalign1.at(3) << endl;
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerRotationZX(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT ZX rotation according to gear file: " << _DUTalign1.at(4) << endl;
                    _DUTalign1.push_back(_siPlanesLayerLayout->getLayerRotationXY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT XY rotation according to gear file: " << _DUTalign1.at(5) << endl;
                }
                if(cPlaneID == _cbcRealDUTsVec.at(1)) {
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerPositionX(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT X shift according to gear file: " << _DUTalign2.at(0) << endl;
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerPositionY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT Y shift according to gear file: " << _DUTalign2.at(1) << endl;
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerPositionZ(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT Z shift according to gear file: " << _DUTalign2.at(2) << endl;
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerRotationZY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT ZY rotation according to gear file: " << _DUTalign2.at(3) << endl;
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerRotationZX(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT ZX rotation according to gear file: " << _DUTalign2.at(4) << endl;
                    _DUTalign2.push_back(_siPlanesLayerLayout->getLayerRotationXY(i));
                    streamlog_out ( MESSAGE5 ) << "Set DUT XY rotation according to gear file: " << _DUTalign2.at(5) << endl;
                }


	}

        // load alignment from file
        _alignmentloaded = false;
        _pre_alignmentloaded = false;
        _referenceloaded = false;

        // allocate arrays for alignment constants
        _dut_align_x_1 = new double[_alignmentCollectionName.size()];
        _dut_align_y_1 = new double[_alignmentCollectionName.size()];
        _dut_align_z_1 = new double[_alignmentCollectionName.size()];
        _dut_align_a_1 = new double[_alignmentCollectionName.size()];
        _dut_align_b_1 = new double[_alignmentCollectionName.size()];
        _dut_align_c_1 = new double[_alignmentCollectionName.size()];
        _dut_align_x_error_1 = new double[_alignmentCollectionName.size()];
        _dut_align_y_error_1 = new double[_alignmentCollectionName.size()];
        _dut_align_z_error_1 = new double[_alignmentCollectionName.size()];
        _dut_align_a_error_1 = new double[_alignmentCollectionName.size()];
        _dut_align_b_error_1 = new double[_alignmentCollectionName.size()];
        _dut_align_c_error_1 = new double[_alignmentCollectionName.size()];
        _dut_pre_align_x_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_y_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_z_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_a_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_b_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_c_1 = new double[_pre_alignmentCollectionName.size()];
        _dut_align_x_2 = new double[_alignmentCollectionName.size()];
        _dut_align_y_2 = new double[_alignmentCollectionName.size()];
        _dut_align_z_2 = new double[_alignmentCollectionName.size()];
        _dut_align_a_2 = new double[_alignmentCollectionName.size()];
        _dut_align_b_2 = new double[_alignmentCollectionName.size()];
        _dut_align_c_2 = new double[_alignmentCollectionName.size()];
        _dut_align_x_error_2 = new double[_alignmentCollectionName.size()];
        _dut_align_y_error_2 = new double[_alignmentCollectionName.size()];
        _dut_align_z_error_2 = new double[_alignmentCollectionName.size()];
        _dut_align_a_error_2 = new double[_alignmentCollectionName.size()];
        _dut_align_b_error_2 = new double[_alignmentCollectionName.size()];
        _dut_align_c_error_2 = new double[_alignmentCollectionName.size()];
        _dut_pre_align_x_2 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_y_2 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_z_2 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_a_2 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_b_2 = new double[_pre_alignmentCollectionName.size()];
        _dut_pre_align_c_2 = new double[_pre_alignmentCollectionName.size()];
}


void CBCHitRecovery::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    bookHistos ( );

}


void CBCHitRecovery::processEvent ( LCEvent * anEvent )
{

    // load the alignment only once
    if (_alignmentloaded == false)
    {
            getAlignment(anEvent);
            _alignmentloaded = true;
            streamlog_out ( MESSAGE5 ) <<  "Alignment loaded..." << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }

    // also load prealignment only once
    if (_pre_alignmentloaded == false)
    {
            getPreAlignment(anEvent);
            _pre_alignmentloaded = true;
            streamlog_out ( MESSAGE5 ) <<  "Prealignment loaded..." << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }

    // also load reference collection only once
    if (_referenceloaded == false)
    {
            getReference(anEvent);
            _referenceloaded = true;
            streamlog_out ( MESSAGE5 ) <<  "Reference loaded..." << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
            streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }


    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

	// the collection we read
	LCCollectionVec * inputTrackVec;
	LCCollectionVec * inputFitPointVec;
	LCCollectionVec * inputHitsVec;
	
	// the collection we output (fits)
	LCCollectionVec * outputHitCollectionVec;

	try
	{
	    outputHitCollectionVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcFitHitOutputCollectionName ) );
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    outputHitCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
	}

	// the collection we output (tracks)
	LCCollectionVec * outputTrackCollectionVec;

	try
	{
	    outputTrackCollectionVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcTracksOutputCollectionName ) );
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    outputTrackCollectionVec = new LCCollectionVec(LCIO::TRACK);
	}

	// prepare an encoder for the hit collection
	CellIDEncoder < TrackerHitImpl > fitHitEncoder ( EUTELESCOPE::HITENCODING, outputHitCollectionVec );

	// find process the tracks
	try
        {
		// give the collection vec its data
		inputTrackVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _InputTrackCollectionName ) );
		inputFitPointVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _InputFitHitsCollectionName ) );
		inputHitsVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcInputCollectionName ) );
		
		// getting the DUT normal vector in the first event
		if( anEvent -> getEventNumber ( ) == 0 ) {
			// get dut normal
			LCCollectionVec* dutnormalvec_collection = dynamic_cast < LCCollectionVec * > (anEvent->getCollection("dutnormal"));
			// get object
			LCGenericObjectImpl* dutnormalvec = dynamic_cast < LCGenericObjectImpl * > (dutnormalvec_collection->at(0));

			// allocare
			_VirtualDutNormal_zplus = new double[3];
			_VirtualDutNormal_zminus = new double[3];
			_VirtualDutPos = new double[3];

			// dut normal vec			
			for(int i = 0; i < 3; i++) {
				_VirtualDutNormal_zplus[i] = dutnormalvec->getDoubleVal(i+3);	
				_VirtualDutNormal_zminus[i] = dutnormalvec->getDoubleVal(i+3);	
				_VirtualDutPos[i] = dutnormalvec->getDoubleVal(i)/1000.0;	
			}
			// we need to get the opposite normal vec for pointing upstream
			_VirtualDutNormal_zminus[2] = -1*_VirtualDutNormal_zminus[2];

			streamlog_out ( MESSAGE4 ) << "Loaded the DUT Pos Vector: X = " << _VirtualDutPos[0] << ", Y = " << _VirtualDutPos[1] << ", Z = " << _VirtualDutPos[2] << endl;
			streamlog_out ( MESSAGE4 ) << "Loaded the DUT Normal Vector: X = " << _VirtualDutNormal_zplus[0] << ", Y = " << _VirtualDutNormal_zplus[1] << ", Z = " << _VirtualDutNormal_zplus[2] << endl;
		
			//delete dutnormalvec;
			//delete dutnormalvec_collection;

			// cheatty but easy...
			ofstream normal_file;
			normal_file.open ("normal_vector.dat", std::ofstream::app);
			normal_file << std::fixed << std::setprecision(8);
			normal_file << anEvent->getRunNumber() << "\t" << _VirtualDutNormal_zplus[0] << "\t" << _VirtualDutNormal_zplus[1] << "\t" << _VirtualDutNormal_zplus[2] << "\t" << "\n";
			normal_file.close();
		}

		
		
		int nTracks = inputTrackVec -> getNumberOfElements ( );
		dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["TrackMultiplicity"] ) -> fill (nTracks);		
		/*		
		// loop over the tracks in the event
                for ( int i = 0; i < nTracks; ++i ) 
                {
		    TrackImpl * track = dynamic_cast < TrackImpl * > ( inputTrackVec -> getElementAt ( i ) );

		    //delete track;
		}*/

		// variables
		const double *cTelescope2Pos = 0;
		const double *cVirtualHitPos = 0;
		const double *cTelescope3Pos = 0;

		// find the fit hit points of the track (needed to calculate the fir hits in the dut
		int nEntries = inputFitPointVec -> getNumberOfElements ( );
                for ( int i = 0; i < nEntries; ++i ) 
                {
		    TrackerHitImpl * hit = dynamic_cast < TrackerHitImpl * > ( inputFitPointVec -> getElementAt ( i ) );

        	    // sensor id decoder
	 	    CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputFitPointVec );
            	    int sensorID = inputCellIDDecoder ( hit ) ["sensorID"];
		    // check that we are on the virtual cbc hit
		    if (sensorID == _cbcVirtualDUTId) {
			cVirtualHitPos = hit->getPosition();
			
			dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosX"] ) -> fill (cVirtualHitPos[0]);
			dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosY"] ) -> fill (cVirtualHitPos[1]);
			dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosZ"] ) -> fill (cVirtualHitPos[2]);
		    } else if (sensorID == 2) {
			cTelescope2Pos = hit->getPosition();
		    } else if (sensorID == 3) {
			cTelescope3Pos = hit->getPosition();
		    }

		    //delete hit;
	        }

		// check that there are tracks in the event
		if (nTracks == 1 && nEntries != 0) {
			// and that we have succeeded to get the track positions
			if (cTelescope2Pos == 0 || cVirtualHitPos == 0 || cTelescope3Pos == 0) {
				streamlog_out ( WARNING ) << "No fit hit for some of the planes in the event " <<  anEvent -> getEventNumber ( ) << endl;
			} else {
				// track from upstream (it pointed from dut to the telescope plane 2, then we will reconstruct sensor 60 hit)
				double cTrackVecUpstream[3];
				for(int i = 0; i < 3; i++) cTrackVecUpstream[i] = cTelescope2Pos[i] - cVirtualHitPos[i];

				// track to downstream
				double cTrackVecDownstream[3];
				for(int i = 0; i < 3; i++) cTrackVecDownstream[i] = cTelescope3Pos[i] - cVirtualHitPos[i];

				// print more for debugging
				/*if (anEvent->getEventNumber() < 50) {
					streamlog_out (DEBUG1) << "vec up: " << cTrackVecUpstream[0] << ", " << cTrackVecUpstream[1] << ", " << cTrackVecUpstream[2] << endl;
					streamlog_out (DEBUG1) << "vec down: " << cTrackVecDownstream[0] << ", " << cTrackVecDownstream[1] << ", " << cTrackVecDownstream[2] << endl;
				}*/

				// cos phi between the vectors
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["cosPhi"] ) -> fill (cos_alpha(cTrackVecUpstream,cTrackVecDownstream));

				// cov matrix
				float fitcov[TRKHITNCOVMATRIX] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

				// fit hit pos in dut 0
				double *fit_hit_pos_dut0 = this->hit_pos(cTrackVecUpstream, _VirtualDutNormal_zminus, cVirtualHitPos, -2.0);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosX_DUT0"] ) -> fill (fit_hit_pos_dut0[0]);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosY_DUT0"] ) -> fill (fit_hit_pos_dut0[1]);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosZ_DUT0"] ) -> fill (fit_hit_pos_dut0[2]);

				// push to the collection
				TrackerHitImpl * fitpoint0 = new TrackerHitImpl;
				fitHitEncoder["sensorID"] =  _cbcRealDUTsVec.at(0);
				fitHitEncoder["properties"] = kFittedHit;
				fitHitEncoder.setCellID ( fitpoint0 );
				fitpoint0 -> setPosition ( fit_hit_pos_dut0 );
				fitpoint0 -> setCovMatrix ( fitcov );
				outputHitCollectionVec -> push_back ( fitpoint0 );


				// fit hit pos in dut 1
				double *fit_hit_pos_dut1 = this->hit_pos(cTrackVecDownstream, _VirtualDutNormal_zplus, cVirtualHitPos, 2.0);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosX_DUT1"] ) -> fill (fit_hit_pos_dut1[0]);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosY_DUT1"] ) -> fill (fit_hit_pos_dut1[1]);
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosZ_DUT1"] ) -> fill (fit_hit_pos_dut1[2]);

				// push to the collection
				TrackerHitImpl * fitpoint1 = new TrackerHitImpl;
				fitHitEncoder["sensorID"] =  _cbcRealDUTsVec.at(1);
				fitHitEncoder["properties"] = kFittedHit;
				fitHitEncoder.setCellID ( fitpoint1 );
				fitpoint1 -> setPosition ( fit_hit_pos_dut1 );
				fitpoint1 -> setCovMatrix ( fitcov );
				outputHitCollectionVec -> push_back ( fitpoint1 );

				// for angle scan
				double roi_x_min = -4.0;
				double roi_x_max = 2.0;
				// for threshold scan
				//double roi_x_min = 1.0;
				//double roi_x_max = 7.0;
				// residual cut (um)
                                double cResidualCut = 45;

				// n input hits (all sensors)
				int nHits = inputHitsVec -> getNumberOfElements ( );
				int cTotalHits0 = 0;
				int cTotalHits1 = 0;

                		for ( int iHit = 0; iHit < nHits; ++iHit ) {
					TrackerHitImpl * hit = dynamic_cast < TrackerHitImpl * > ( inputHitsVec -> getElementAt ( iHit ) );

					// sensor id decoder
	 	    			CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputFitPointVec );
            	    			int sensorID = inputCellIDDecoder ( hit ) ["sensorID"];

					// first sensor
					if (sensorID == _cbcRealDUTsVec.at(0)) {
						cTotalHits0++;
					}

					// second sensor
					if (sensorID == _cbcRealDUTsVec.at(1)) {
						cTotalHits1++;
					}
				}

				// tdc of the event 
				int cEventTDC = stoi(anEvent->getParameters().getStringVal("TDC"));


				if(fit_hit_pos_dut0[0] >= roi_x_min && fit_hit_pos_dut0[0] <= roi_x_max) {
					if (cTotalHits0 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["TrackTDC0"] ) -> fill (cEventTDC);
				}
				if(fit_hit_pos_dut1[0] >= roi_x_min && fit_hit_pos_dut1[0] <= roi_x_max) {
					if (cTotalHits1 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["TrackTDC1"] ) -> fill (cEventTDC);
				}

                                streamlog_out ( DEBUG6 ) << endl << endl << _cbcRealDUTsVec.at(0) << ": fit X position: " << fit_hit_pos_dut0[0] << endl;
                                TVector2 hit0 = DoDeAlignment1(TVector3(fit_hit_pos_dut0[0],fit_hit_pos_dut0[1], fit_hit_pos_dut0[2]));
                                if(cTotalHits0 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["TracksPerChannel0"] ) -> fill (roundValue(hit0.X()));

                                streamlog_out ( DEBUG6 ) << endl << endl << _cbcRealDUTsVec.at(1) << ": fit X position: " << fit_hit_pos_dut1[0] << endl;
                                TVector2 hit1 = DoDeAlignment2(TVector3(fit_hit_pos_dut1[0],fit_hit_pos_dut1[1], fit_hit_pos_dut1[2]));
                                if(cTotalHits1 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["TracksPerChannel1"] ) -> fill (roundValue(hit1.X()));


				// now calculate the residuals
                		for ( int iHit = 0; iHit < nHits; ++iHit ) {
					TrackerHitImpl * hit = dynamic_cast < TrackerHitImpl * > ( inputHitsVec -> getElementAt ( iHit ) );

					// sensor id decoder
	 	    			CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputFitPointVec );
            	    			int sensorID = inputCellIDDecoder ( hit ) ["sensorID"];

					// virtual
					if (sensorID == _cbcVirtualDUTId) {
						const double *pos = hit->getPosition();
						double resX = (cVirtualHitPos[0] - pos[0])*1000;
						double resY = (cVirtualHitPos[1] - pos[1])*1000;
						double resZ = (cVirtualHitPos[2] - pos[2])*1000;

						dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualResidualX"] ) -> fill (resX);
					}

					// first sensor
					if (sensorID == _cbcRealDUTsVec.at(0)) {
						const double *pos = hit->getPosition();
                                                TrackerDataImpl* raw_hit = static_cast < TrackerDataImpl* > ( hit -> getRawHits ( )[0]);
                                                EUTelSimpleVirtualCluster * cluster1 = 0;
                                                cluster1 = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( raw_hit );
                                                float x1 = -1.0;
                                                float y1 = -1.0;
                                                if ( cluster1 != 0 )
                                                {
                                                    cluster1 -> getCenterOfGravity ( x1, y1 );
                                                }

                                                double resX = (fit_hit_pos_dut0[0] - pos[0])*1000;
                                                /*double resY = (fit_hit_pos_dut0[1] - pos[1])*1000;
                                                double resZ = (fit_hit_pos_dut0[2] - pos[2])*1000;*/
                                                double resX_dealigned = (hit0.X()-x1)*_pitchx*1000;
                                                if((roundValue(hit0.X())==x1) && cTotalHits0 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["HitsPerChannel0"] ) -> fill (x1);


						if((pos[0] >= roi_x_min && pos[0] <= roi_x_max) && fabs(resX) < cResidualCut) {
							if (cTotalHits0 <= 1) {
								dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["HitTDC0"] ) -> fill (cEventTDC);
							}
						}

                                                streamlog_out ( DEBUG6 ) << _cbcRealDUTsVec.at(0) << ": Global DUT Hit position: (" << pos[0] << " | " << pos[1] << " )" << endl;
                                                streamlog_out ( DEBUG6 ) << _cbcRealDUTsVec.at(0) << ": Local DUT Hit position: (" << x1 << " | " << y1 << " )" << endl;


                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT0"] ) -> fill (resX);
                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT0_DeAligned"] ) -> fill (resX_dealigned);
                                        }

					// second sensor
					if (sensorID == _cbcRealDUTsVec.at(1)) {
						const double *pos = hit->getPosition();
                                                TrackerDataImpl* raw_hit = static_cast < TrackerDataImpl* > ( hit -> getRawHits ( )[0]);
                                                EUTelSimpleVirtualCluster * cluster1 = 0;
                                                cluster1 = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( raw_hit );
                                                float x1 = -1.0;
                                                float y1 = -1.0;
                                                if ( cluster1 != 0 )
                                                {
                                                    cluster1 -> getCenterOfGravity ( x1, y1 );
                                                }

                                                double resX = (fit_hit_pos_dut1[0] - pos[0])*1000;
                                                /*double resY = (fit_hit_pos_dut1[1] - pos[1])*1000;
                                                double resZ = (fit_hit_pos_dut1[2] - pos[2])*1000;*/
                                                double resX_dealigned = (hit1.X()-x1)*_pitchx*1000;
                                                if((roundValue(hit1.X())==x1) && cTotalHits1 <= 1) dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["HitsPerChannel1"] ) -> fill (x1);


						if((pos[0] >= roi_x_min && pos[0] <= roi_x_max) && fabs(resX) < cResidualCut) {
							if (cTotalHits1 <= 1) {
								dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["HitTDC1"] ) -> fill (cEventTDC);
							}
                                                }

                                                streamlog_out ( DEBUG6 ) << _cbcRealDUTsVec.at(1) << ": Global DUT Hit position: (" << pos[0] << " | " << pos[1] << " )" << endl;
                                                streamlog_out ( DEBUG6 ) << _cbcRealDUTsVec.at(1) << ": Local DUT Hit position: (" << x1 << " | " << y1 << " )" << endl;

                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT1"] ) -> fill (resX);
                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT1_DeAligned"] ) -> fill (resX_dealigned);

					}
				}

				// delete the hist positions
				delete fit_hit_pos_dut0;
				delete fit_hit_pos_dut1;

				// clear pos vectors
				//delete cTelescope2Pos;
				//delete cVirtualHitPos;
				//delete cTelescope3Pos;
			}
		}

	}
	catch ( lcio::DataNotAvailableException )
	{

	}

	anEvent->addCollection( outputHitCollectionVec, _cbcFitHitOutputCollectionName );
	anEvent->addCollection( outputTrackCollectionVec, _cbcTracksOutputCollectionName );

}


void CBCHitRecovery::check ( LCEvent * /* evt */ )
{

}


void CBCHitRecovery::end ( )
{
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;

}


void CBCHitRecovery::fillHistos ( )
{

}


void CBCHitRecovery::bookHistos ( )
{
	string basePath = "Input";
	AIDAProcessor::tree ( this ) -> mkdir ( basePath.c_str ( ) );
	basePath.append ( "/" );

	AIDA::IHistogram1D * cTrackMultiplicity = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "TrackMultiplicity").c_str(), 15, 0, 15 );
	_aidaHistoMap.insert ( make_pair ( "TrackMultiplicity", cTrackMultiplicity ) );
	cTrackMultiplicity -> setTitle ( "Number of tracks per Event;N tracks;Entries" );

	AIDA::IHistogram1D * cVirtualHitXHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosX" ).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "VirtualHitPosX", cVirtualHitXHist ) );
	cVirtualHitXHist -> setTitle ( "Virtual Hit Position X;X [mm];Entries" );

	AIDA::IHistogram1D * cVirtualHitYHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosY" ).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "VirtualHitPosY", cVirtualHitYHist ) );
	cVirtualHitYHist -> setTitle ( "Virtual Hit Position Y;Y [mm];Entries" );

	AIDA::IHistogram1D * cVirtualHitZHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosZ" ).c_str ( ), 100, 365, 375 );
	_aidaHistoMap.insert ( make_pair ( "VirtualHitPosZ", cVirtualHitZHist ) );
	cVirtualHitZHist -> setTitle ( "Virtual Hit Position Z;Z [mm];Entries" );

	AIDA::IHistogram1D * cCosPhiHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "cosPhi" ).c_str ( ), 500, -1.000001, -0.99999 );
	_aidaHistoMap.insert ( make_pair ( "cosPhi", cCosPhiHist ) );
	cCosPhiHist -> setTitle ( "cos between inbound and outbound track vectors;cos(phi);Entries" );

	string basePathOut = "Output";
	AIDAProcessor::tree ( this ) -> mkdir ( basePathOut.c_str ( ) );
	basePathOut.append ( "/" );

	AIDA::IHistogram1D * cDUT0FitHitPosX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosX_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosX_DUT0", cDUT0FitHitPosX ) );
	cDUT0FitHitPosX -> setTitle ( "Fit Hit Pos X - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";X [mm];Entries" );

	AIDA::IHistogram1D * cDUT0FitHitPosY = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosY_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosY_DUT0", cDUT0FitHitPosY ) );
	cDUT0FitHitPosY -> setTitle ( "Fit Hit Pos Y - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";Y [mm];Entries" );

	AIDA::IHistogram1D * cDUT0FitHitPosZ = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosZ_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, 365, 375 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosZ_DUT0", cDUT0FitHitPosZ ) );
	cDUT0FitHitPosX -> setTitle ( "Fit Hit Pos Z - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";Z [mm];Entries" );

	AIDA::IHistogram1D * cDUT1FitHitPosX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosX_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosX_DUT1", cDUT1FitHitPosX ) );
	cDUT1FitHitPosX -> setTitle ( "Fit Hit Pos X - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";X [mm];Entries" );

	AIDA::IHistogram1D * cDUT1FitHitPosY = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosY_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, -100, 100 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosY_DUT1", cDUT1FitHitPosY ) );
	cDUT1FitHitPosY -> setTitle ( "Fit Hit Pos Y - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";Y [mm];Entries" );

	AIDA::IHistogram1D * cDUT1FitHitPosZ = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosZ_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, 365, 375 );
	_aidaHistoMap.insert ( make_pair ( "FitHitPosZ_DUT1", cDUT1FitHitPosZ ) );
	cDUT1FitHitPosX -> setTitle ( "Fit Hit Pos Z - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";Z [mm];Entries" );

	AIDA::IHistogram1D * cVirtualResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "VirtualResidualX" ).c_str ( ), 100, -500, 500 );
	_aidaHistoMap.insert ( make_pair ( "VirtualResidualX", cVirtualResidualX ) );
	cVirtualResidualX -> setTitle ( "Residual X - Virtual DUT;X [um];Entries" );

        AIDA::IHistogram1D * cDUT0ResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 600, -1500, 1500 );
	_aidaHistoMap.insert ( make_pair ( "ResidualX_DUT0", cDUT0ResidualX ) );
	cDUT0ResidualX -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";X [um];Entries" );

        AIDA::IHistogram1D * cDUT1ResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 600, -1500, 1500 );
	_aidaHistoMap.insert ( make_pair ( "ResidualX_DUT1", cDUT1ResidualX ) );
	cDUT1ResidualX -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";X [um];Entries" );

        AIDA::IHistogram1D * cDUT0ResidualX_DeAligned = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(0)) + "_DeAligned").c_str ( ), 600, -1500, 1500 );
        _aidaHistoMap.insert ( make_pair ( "ResidualX_DUT0_DeAligned", cDUT0ResidualX_DeAligned ) );
        cDUT0ResidualX_DeAligned -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";X [um];Entries" );

        AIDA::IHistogram1D * cDUT1ResidualX_DeAligned = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(1)) + "_DeAligned").c_str ( ), 600, -1500, 1500 );
        _aidaHistoMap.insert ( make_pair ( "ResidualX_DUT1_DeAligned", cDUT1ResidualX_DeAligned ) );
        cDUT1ResidualX_DeAligned -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";X [um];Entries" );

	AIDA::IHistogram1D * cTDCHit0 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "HitTDC" + to_string(_cbcRealDUTsVec.at(0))).c_str(), 8, 0, 8 );
	_aidaHistoMap.insert ( make_pair ( "HitTDC0", cTDCHit0 ) );
	cTDCHit0 -> setTitle ( "Hit TDC;TDC;Entries" );

	AIDA::IHistogram1D * cTDCHit1 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "HitTDC" + to_string(_cbcRealDUTsVec.at(1))).c_str(), 8, 0, 8 );
	_aidaHistoMap.insert ( make_pair ( "HitTDC1", cTDCHit1 ) );
	cTDCHit1 -> setTitle ( "Hit TDC;TDC;Entries" );

	AIDA::IHistogram1D * cTDCTrack0 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "TrackTDC" + to_string(_cbcRealDUTsVec.at(0))).c_str(), 8, 0, 8 );
	_aidaHistoMap.insert ( make_pair ( "TrackTDC0", cTDCTrack0 ) );
	cTDCTrack0 -> setTitle ( "Track TDC;TDC;Entries" );

	AIDA::IHistogram1D * cTDCTrack1 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "TrackTDC" + to_string(_cbcRealDUTsVec.at(1))).c_str(), 8, 0, 8 );
	_aidaHistoMap.insert ( make_pair ( "TrackTDC1", cTDCTrack1 ) );
	cTDCTrack1 -> setTitle ( "Track TDC;TDC;Entries" );

        AIDA::IHistogram1D * cTracksPerCh0 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "TracksPerChannel" + to_string(_cbcRealDUTsVec.at(0))).c_str(), 254, 0, 254 );
        _aidaHistoMap.insert ( make_pair ( "TracksPerChannel0", cTracksPerCh0 ) );
        cTracksPerCh0 -> setTitle ( "Tracks Per Channel;Channel;Entries" );

        AIDA::IHistogram1D * cTracksPerCh1 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "TracksPerChannel" + to_string(_cbcRealDUTsVec.at(1))).c_str(), 254, 0, 254 );
        _aidaHistoMap.insert ( make_pair ( "TracksPerChannel1", cTracksPerCh1 ) );
        cTracksPerCh1 -> setTitle ( "Tracks Per Channel;Channel;Entries" );

        AIDA::IHistogram1D * cHitsPerCh0 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "HitsPerChannel" + to_string(_cbcRealDUTsVec.at(0))).c_str(), 254, 0, 254 );
        _aidaHistoMap.insert ( make_pair ( "HitsPerChannel0", cHitsPerCh0 ) );
        cHitsPerCh0 -> setTitle ( "Hits Per Channel;Channel;Entries" );

        AIDA::IHistogram1D * cHitsPerCh1 = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "HitsPerChannel" + to_string(_cbcRealDUTsVec.at(1))).c_str(), 254, 0, 254 );
        _aidaHistoMap.insert ( make_pair ( "HitsPerChannel1", cHitsPerCh1 ) );
        cHitsPerCh1 -> setTitle ( "Hits Per Channel;Channel;Entries" );

}

// get the module of vector
double CBCHitRecovery::vector_get_length(const double *vec) {
	double length = 0;
	for(int i = 0; i < 3; i++) length += vec[i]*vec[i];
	return sqrt(length);
}

// set certain length of the vector
void CBCHitRecovery::vector_set_length(double *& vec, double length) {
	// we also need to normalize it
	double current_length = vector_get_length(vec);
	// now set the new length
	for(int i = 0; i < 3; i++) vec[i] = vec[i]*length/current_length;
}

// angle between two vector
double CBCHitRecovery::cos_alpha(const double *vec1, const double *vec2) {
	// variable 
	double angle = 0;
	// calculate the scalar
	for(int i = 0; i < 3; i++) angle += vec1[i]*vec2[i];
	// now divide by modules
	angle = angle/(vector_get_length(vec1)*vector_get_length(vec2));
	// return
	return angle;	
}

// finds the position of fithit in real plane
double* CBCHitRecovery::hit_pos(const double track[],const double normal[],const double virtual_hit[], double distance) {
	// init vector
	double *pos = new double[3];
	// module of distance
	double abs_distance = distance;
	if(abs_distance < 0) abs_distance = -abs_distance;
	// create copies of the vectors to modify them
	double *track1 = new double[3];
	for(int i = 0; i < 3; i++) {
		track1[i] = track[i];
	}
	// set the length of the track vector to distance/cos(alpha), alpha - angle between normal and track
	vector_set_length(track1, abs_distance/cos_alpha(track, normal));
	// add two to get the position
	for(int i = 0; i < 3; i++) pos[i] = track1[i] + virtual_hit[i];
	// clear
	delete track1;
	// return
	return pos;
}

TVector2 CBCHitRecovery::DoDeAlignment1(TVector3 hit_input) {

    // store the global position in a root 3vectot (get a copy)
    TVector3 iterationvector = hit_input;

    // dealign the reference hit too
    double refx = 0.0;
    double refy = 0.0;
    double refz = 0.0;

    // the central coordinates of our DUT:
    double dutcenter_x = 0.0;
    double dutcenter_y = 0.0;
    double dutcenter_z = 0.0;

    // undo the alignment
    for (int i = _alignmentCollectionName.size(); i>0 ; i--)
    {
            streamlog_out ( DEBUG2 ) << "Undoing alignment, step " << i << " !" << endl;

            // i-1 since the array starts at 0!

            // align the reference
            refx = _x_refhit_1 + _dut_align_x_1[i-1];
            refy = _y_refhit_1 + _dut_align_y_1[i-1];
            refz = _z_refhit_1 + _dut_align_z_1[i-1];

            TVector3 _RotatedVector( _a_refhit_1, _b_refhit_1, _c_refhit_1 );
            _RotatedVector.RotateZ( _dut_align_a_1[i-1] );
            _RotatedVector.RotateY( _dut_align_b_1[i-1] );
            _RotatedVector.RotateX( _dut_align_c_1[i-1] );

            // align our global track (iterationvector)
            // relative to the reference
            double inputPosition[3] = { iterationvector.X(), iterationvector.Y(), iterationvector.Z() };

            dutcenter_x = refx;
            dutcenter_y = refy;
            dutcenter_z = refz;

            dutcenter_x -= _dut_align_x_1[i-1];
            dutcenter_y -= _dut_align_y_1[i-1];
            dutcenter_z -= _dut_align_z_1[i-1];

            inputPosition[0] = inputPosition[0] - dutcenter_x;
            inputPosition[1] = inputPosition[1] - dutcenter_y;
            inputPosition[2] = inputPosition[2] - dutcenter_z;

            double outputPosition[3] = {0.0, 0.0, 0.0};

            outputPosition[0] = dutcenter_x;
            outputPosition[1] = dutcenter_y;
            outputPosition[2] = dutcenter_z;

            TVector3 iCenterOfSensorFrame(  inputPosition[0],  inputPosition[1],  inputPosition[2] );

            iCenterOfSensorFrame.RotateZ( _dut_align_c_1[i-1] );
            iCenterOfSensorFrame.RotateY( _dut_align_b_1[i-1] );
            iCenterOfSensorFrame.RotateX( _dut_align_a_1[i-1] );

            outputPosition[0] += iCenterOfSensorFrame(0);
            outputPosition[1] += iCenterOfSensorFrame(1);
            outputPosition[2] += iCenterOfSensorFrame(2);

            // second the shift
            outputPosition[0] += _dut_align_x_1[i-1];
            outputPosition[1] += _dut_align_y_1[i-1];
            outputPosition[2] += _dut_align_z_1[i-1];

            iterationvector.SetX(outputPosition[0]);
            iterationvector.SetY(outputPosition[1]);
            iterationvector.SetZ(outputPosition[2]);


            streamlog_out ( DEBUG2 ) << "De-aligned vector is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // undo the prealignment
    for (int i=_pre_alignmentCollectionName.size();i>0;i--)
    {
            streamlog_out ( DEBUG2 ) << "Undoing prealignment, step " << i << " !" << endl;

            // i-1 since the array starts at 0!

            // align the reference
            refx = _x_refhit_1 + _dut_pre_align_x_1[i-1];
            refy = _y_refhit_1 + _dut_pre_align_y_1[i-1];
            refz = _z_refhit_1 + _dut_pre_align_z_1[i-1];

            TVector3 _RotatedVector( _a_refhit_1, _b_refhit_1, _c_refhit_1 );
            _RotatedVector.RotateZ( _dut_pre_align_a_1[i-1] );
            _RotatedVector.RotateY( _dut_pre_align_b_1[i-1] );
            _RotatedVector.RotateX( _dut_pre_align_c_1[i-1] );

            // align our global track (iterationvector)
            // relative to the reference
            double inputPosition[3] = { iterationvector.X(), iterationvector.Y(), iterationvector.Z() };

            dutcenter_x = refx;
            dutcenter_y = refy;
            dutcenter_z = refz;

            dutcenter_x -= _dut_pre_align_x_1[i-1];
            dutcenter_y -= _dut_pre_align_y_1[i-1];
            dutcenter_z -= _dut_pre_align_z_1[i-1];

            inputPosition[0] = inputPosition[0] - dutcenter_x;
            inputPosition[1] = inputPosition[1] - dutcenter_y;
            inputPosition[2] = inputPosition[2] - dutcenter_z;

            double outputPosition[3] = {0.0, 0.0, 0.0};

            outputPosition[0] = dutcenter_x;
            outputPosition[1] = dutcenter_y;
            outputPosition[2] = dutcenter_z;

            TVector3 iCenterOfSensorFrame(  inputPosition[0],  inputPosition[1],  inputPosition[2] );

            iCenterOfSensorFrame.RotateZ( _dut_pre_align_c_1[i-1] );
            iCenterOfSensorFrame.RotateY( _dut_pre_align_b_1[i-1] );
            iCenterOfSensorFrame.RotateX( _dut_pre_align_a_1[i-1] );

            outputPosition[0] += iCenterOfSensorFrame(0);
            outputPosition[1] += iCenterOfSensorFrame(1);
            outputPosition[2] += iCenterOfSensorFrame(2);

            // second the shift
            outputPosition[0] += _dut_pre_align_x_1[i-1];
            outputPosition[1] += _dut_pre_align_y_1[i-1];
            outputPosition[2] += _dut_pre_align_z_1[i-1];

            iterationvector.SetX(outputPosition[0]);
            iterationvector.SetY(outputPosition[1]);
            iterationvector.SetZ(outputPosition[2]);

            streamlog_out ( DEBUG2 ) << "De-prealigned vector is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // FIXME only x and y...
    if (_useOriginalPreAlignment)
    {
            iterationvector.SetX(iterationvector.X() + _dut_original_pre_align_x_1);
            iterationvector.SetY(iterationvector.Y() + _dut_original_pre_align_y_1);

            streamlog_out ( DEBUG2 ) << "Vector after undoing original prealignment is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // iterationvector is now in local coordinates, still in mm
    // revert gear rotations and shifts too
    // FIXME for now only alpha

    float tempx = iterationvector.X();
    float rotx = tempx/cos(_DUTalign1.at(4)*PI/180.0);
    iterationvector.SetX(rotx-_DUTalign1.at(0));

    streamlog_out ( DEBUG6 ) << " " << endl;
    streamlog_out ( DEBUG6 ) << "Local track hit:  ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z() << " )" << endl;

    double dutTrackX_local_pix = (_pixelx / 2.0) - 0.5 - iterationvector.X() / _pitchx;
    double dutTrackY_local_pix = (_pixely / 2.0) - 0.5 - iterationvector.Y() / _pitchy;

    streamlog_out ( DEBUG6 ) << "Points to pixels: ( " << dutTrackX_local_pix << " | " << dutTrackY_local_pix << " )" << endl;
    streamlog_out ( DEBUG6 ) << " " << endl;

    return TVector2(dutTrackX_local_pix,dutTrackY_local_pix);

}

TVector2 CBCHitRecovery::DoDeAlignment2(TVector3 hit_input) {

    // store the global position in a root 3vectot (get a copy)
    TVector3 iterationvector = hit_input;

    // dealign the reference hit too
    double refx = 0.0;
    double refy = 0.0;
    double refz = 0.0;

    // the central coordinates of our DUT:
    double dutcenter_x = 0.0;
    double dutcenter_y = 0.0;
    double dutcenter_z = 0.0;

    // undo the alignment
    for (int i = _alignmentCollectionName.size(); i>0 ; i--)
    {
            streamlog_out ( DEBUG2 ) << "Undoing alignment, step " << i << " !" << endl;

            // i-1 since the array starts at 0!

            // align the reference
            refx = _x_refhit_2 + _dut_align_x_2[i-1];
            refy = _y_refhit_2 + _dut_align_y_2[i-1];
            refz = _z_refhit_2 + _dut_align_z_2[i-1];

            TVector3 _RotatedVector( _a_refhit_2, _b_refhit_2, _c_refhit_2 );
            _RotatedVector.RotateZ( _dut_align_a_2[i-1] );
            _RotatedVector.RotateY( _dut_align_b_2[i-1] );
            _RotatedVector.RotateX( _dut_align_c_2[i-1] );

            // align our global track (iterationvector)
            // relative to the reference
            double inputPosition[3] = { iterationvector.X(), iterationvector.Y(), iterationvector.Z() };

            dutcenter_x = refx;
            dutcenter_y = refy;
            dutcenter_z = refz;

            dutcenter_x -= _dut_align_x_2[i-1];
            dutcenter_y -= _dut_align_y_2[i-1];
            dutcenter_z -= _dut_align_z_2[i-1];

            inputPosition[0] = inputPosition[0] - dutcenter_x;
            inputPosition[1] = inputPosition[1] - dutcenter_y;
            inputPosition[2] = inputPosition[2] - dutcenter_z;

            double outputPosition[3] = {0.0, 0.0, 0.0};

            outputPosition[0] = dutcenter_x;
            outputPosition[1] = dutcenter_y;
            outputPosition[2] = dutcenter_z;

            TVector3 iCenterOfSensorFrame(  inputPosition[0],  inputPosition[1],  inputPosition[2] );

            iCenterOfSensorFrame.RotateZ( _dut_align_c_2[i-1] );
            iCenterOfSensorFrame.RotateY( _dut_align_b_2[i-1] );
            iCenterOfSensorFrame.RotateX( _dut_align_a_2[i-1] );

            outputPosition[0] += iCenterOfSensorFrame(0);
            outputPosition[1] += iCenterOfSensorFrame(1);
            outputPosition[2] += iCenterOfSensorFrame(2);

            // second the shift
            outputPosition[0] += _dut_align_x_2[i-1];
            outputPosition[1] += _dut_align_y_2[i-1];
            outputPosition[2] += _dut_align_z_2[i-1];

            iterationvector.SetX(outputPosition[0]);
            iterationvector.SetY(outputPosition[1]);
            iterationvector.SetZ(outputPosition[2]);


            streamlog_out ( DEBUG2 ) << "De-aligned vector is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // undo the prealignment
    for (int i=_pre_alignmentCollectionName.size();i>0;i--)
    {
            streamlog_out ( DEBUG2 ) << "Undoing prealignment, step " << i << " !" << endl;

            // i-1 since the array starts at 0!

            // align the reference
            refx = _x_refhit_2 + _dut_pre_align_x_2[i-1];
            refy = _y_refhit_2 + _dut_pre_align_y_2[i-1];
            refz = _z_refhit_2 + _dut_pre_align_z_2[i-1];

            TVector3 _RotatedVector( _a_refhit_2, _b_refhit_2, _c_refhit_2 );
            _RotatedVector.RotateZ( _dut_pre_align_a_2[i-1] );
            _RotatedVector.RotateY( _dut_pre_align_b_2[i-1] );
            _RotatedVector.RotateX( _dut_pre_align_c_2[i-1] );

            // align our global track (iterationvector)
            // relative to the reference
            double inputPosition[3] = { iterationvector.X(), iterationvector.Y(), iterationvector.Z() };

            dutcenter_x = refx;
            dutcenter_y = refy;
            dutcenter_z = refz;

            dutcenter_x -= _dut_pre_align_x_2[i-1];
            dutcenter_y -= _dut_pre_align_y_2[i-1];
            dutcenter_z -= _dut_pre_align_z_2[i-1];

            inputPosition[0] = inputPosition[0] - dutcenter_x;
            inputPosition[1] = inputPosition[1] - dutcenter_y;
            inputPosition[2] = inputPosition[2] - dutcenter_z;

            double outputPosition[3] = {0.0, 0.0, 0.0};

            outputPosition[0] = dutcenter_x;
            outputPosition[1] = dutcenter_y;
            outputPosition[2] = dutcenter_z;

            TVector3 iCenterOfSensorFrame(  inputPosition[0],  inputPosition[1],  inputPosition[2] );

            iCenterOfSensorFrame.RotateZ( _dut_pre_align_c_2[i-1] );
            iCenterOfSensorFrame.RotateY( _dut_pre_align_b_2[i-1] );
            iCenterOfSensorFrame.RotateX( _dut_pre_align_a_2[i-1] );

            outputPosition[0] += iCenterOfSensorFrame(0);
            outputPosition[1] += iCenterOfSensorFrame(1);
            outputPosition[2] += iCenterOfSensorFrame(2);

            // second the shift
            outputPosition[0] += _dut_pre_align_x_2[i-1];
            outputPosition[1] += _dut_pre_align_y_2[i-1];
            outputPosition[2] += _dut_pre_align_z_2[i-1];

            iterationvector.SetX(outputPosition[0]);
            iterationvector.SetY(outputPosition[1]);
            iterationvector.SetZ(outputPosition[2]);

            streamlog_out ( DEBUG2 ) << "De-prealigned vector is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // FIXME only x and y...
    if (_useOriginalPreAlignment)
    {
            iterationvector.SetX(iterationvector.X() + _dut_original_pre_align_x_2);
            iterationvector.SetY(iterationvector.Y() + _dut_original_pre_align_y_2);

            streamlog_out ( DEBUG2 ) << "Vector after undoing original prealignment is:   ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z()<< " ) !" << endl;
    }

    // iterationvector is now in local coordinates, still in mm
    // revert gear rotations and shifts too
    // FIXME for now only alpha

    float tempx = iterationvector.X();
    float rotx = tempx/cos(_DUTalign2.at(4)*PI/180.0);
    iterationvector.SetX(rotx-_DUTalign2.at(0));

    streamlog_out ( DEBUG6 ) << " " << endl;
    streamlog_out ( DEBUG6 ) << "Local track hit:  ( " << iterationvector.X() << " | " << iterationvector.Y() << " | " << iterationvector.Z() << " )" << endl;

    double dutTrackX_local_pix = (_pixelx / 2.0) - 0.5 - iterationvector.X() / _pitchx;
    double dutTrackY_local_pix = (_pixely / 2.0) - 0.5 - iterationvector.Y() / _pitchy;

    streamlog_out ( DEBUG6 ) << "Points to pixels: ( " << dutTrackX_local_pix << " | " << dutTrackY_local_pix << " )" << endl;
    streamlog_out ( DEBUG6 ) << " " << endl;

    return TVector2(dutTrackX_local_pix,dutTrackY_local_pix);

}

void CBCHitRecovery::getReference (LCEvent * event)
{

        // get the DUT z position from the reference file
        // there might be a difference between this z and the one loaded from the gear file

        _x_refhit_1 = 0.0;
        _y_refhit_1 = 0.0;
        _z_refhit_1 = 0.0;
        _a_refhit_1 = 0.0;
        _b_refhit_1 = 0.0;
        _c_refhit_1 = 0.0;

        _x_refhit_2 = 0.0;
        _y_refhit_2 = 0.0;
        _z_refhit_2 = 0.0;
        _a_refhit_2 = 0.0;
        _b_refhit_2 = 0.0;
        _c_refhit_2 = 0.0;

        LCCollectionVec* referenceHitVec = dynamic_cast < LCCollectionVec * > ( event->getCollection( _referencecollectionname) );
        for(size_t i = 0 ; i <  static_cast< size_t >(referenceHitVec->getNumberOfElements()); i++)
        {
                streamlog_out ( DEBUG7 ) << "Looping references... " << i << endl;
                EUTelReferenceHit * refhit = static_cast< EUTelReferenceHit*> ( referenceHitVec->getElementAt(i) );
                if( _cbcRealDUTsVec.at(0) == refhit->getSensorID() )
                {
                        _x_refhit_1 = refhit->getXOffset();
                        _y_refhit_1 = refhit->getYOffset();
                        _z_refhit_1 = refhit->getZOffset();
                        _a_refhit_1 = refhit->getAlpha();
                        _b_refhit_1 = refhit->getBeta();
                        _c_refhit_1 = refhit->getGamma();
                }

                if( _cbcRealDUTsVec.at(1) == refhit->getSensorID() )
                {
                        _x_refhit_2 = refhit->getXOffset();
                        _y_refhit_2 = refhit->getYOffset();
                        _z_refhit_2 = refhit->getZOffset();
                        _a_refhit_2 = refhit->getAlpha();
                        _b_refhit_2 = refhit->getBeta();
                        _c_refhit_2 = refhit->getGamma();
                }
        }

        streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
        streamlog_out ( DEBUG7 ) << "Reference Collection loaded:" << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference X is:     " << _x_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference Y is:     " << _y_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference Z is:     " << _z_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference Alpha is: " << _a_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference Beta is:  " << _b_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT1 reference Gamma is: " << _c_refhit_1 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference X is:     " << _x_refhit_2 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference Y is:     " << _y_refhit_2 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference Z is:     " << _z_refhit_2 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference Alpha is: " << _a_refhit_2 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference Beta is:  " << _b_refhit_2 << endl;
        streamlog_out ( DEBUG7 ) << "DUT2 reference Gamma is: " << _c_refhit_2 << endl;


}

void CBCHitRecovery::getAlignment (LCEvent * event)
{

        // Get the alignment from file
        // assumes rotations in rad!

        LCCollectionVec * alignmentCollection;

        for (size_t i =0; i<_alignmentCollectionName.size(); i++)
        {

                try {
                        alignmentCollection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _alignmentCollectionName[i] ) );

                        for (int j=0; j<_nTelPlanes;j++)
                        {
                                 EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollection->getElementAt( j ) );

                                 // fixme has to be _cbcRealDUTsVec
                                 if (alignment->getSensorID() == _cbcVirtualDUTId)
                                 {
                                        _dut_align_x_1[i] = alignment->getXOffset();
                                        _dut_align_y_1[i] = alignment->getYOffset();
                                        _dut_align_z_1[i] = alignment->getZOffset();
                                        _dut_align_a_1[i] = alignment->getAlpha();
                                        _dut_align_b_1[i] = alignment->getBeta();
                                        _dut_align_c_1[i] = alignment->getGamma();
                                        _dut_align_x_error_1[i] = alignment->getXOffsetError();
                                        _dut_align_y_error_1[i] = alignment->getYOffsetError();
                                        _dut_align_z_error_1[i] = alignment->getZOffsetError();
                                        _dut_align_a_error_1[i] = alignment->getAlphaError();
                                        _dut_align_b_error_1[i] = alignment->getBetaError();
                                        _dut_align_c_error_1[i] = alignment->getGammaError();
                                }
                                 // fixme has to be _cbcRealDUTsVec
                                if (alignment->getSensorID() == _cbcVirtualDUTId)
                                {
                                    _dut_align_x_2[i] = alignment->getXOffset();
                                    _dut_align_y_2[i] = alignment->getYOffset();
                                    _dut_align_z_2[i] = alignment->getZOffset();
                                    _dut_align_a_2[i] = alignment->getAlpha();
                                    _dut_align_b_2[i] = alignment->getBeta();
                                    _dut_align_c_2[i] = alignment->getGamma();
                                    _dut_align_x_error_2[i] = alignment->getXOffsetError();
                                    _dut_align_y_error_2[i] = alignment->getYOffsetError();
                                    _dut_align_z_error_2[i] = alignment->getZOffsetError();
                                    _dut_align_a_error_2[i] = alignment->getAlphaError();
                                    _dut_align_b_error_2[i] = alignment->getBetaError();
                                    _dut_align_c_error_2[i] = alignment->getGammaError();
                                }
                        }
                } catch (lcio::DataNotAvailableException& e) {
                        streamlog_out( ERROR5 ) << "No alignment collection with name " << _alignmentCollectionName[i] << " !" << endl;
                }

                streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
                streamlog_out ( DEBUG7 ) << "Alignment loading, iteration " << i+1 << ":" << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: X is:     " << _dut_align_x_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: Y is:     " << _dut_align_y_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: Z is:     " << _dut_align_z_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: Alpha is: " << _dut_align_a_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: Beta is:  " << _dut_align_b_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 alignment: Gamma is: " << _dut_align_c_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: X is:     " << _dut_align_x_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: Y is:     " << _dut_align_y_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: Z is:     " << _dut_align_z_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: Alpha is: " << _dut_align_a_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: Beta is:  " << _dut_align_b_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 alignment: Gamma is: " << _dut_align_c_2[i] << endl;
        }
}

void CBCHitRecovery::getPreAlignment (LCEvent * event)
{

        // Get the alignment from file
        // assumes rotations in rad!

        LCCollectionVec * pre_alignmentCollection;

        for (size_t i =0; i<_pre_alignmentCollectionName.size(); i++)
        {
                try {
                        pre_alignmentCollection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _pre_alignmentCollectionName[i] ) );

                        for (int j=0; j<_nTelPlanes;j++)
                        {
                                 EUTelAlignmentConstant * pre_alignment = static_cast< EUTelAlignmentConstant * > ( pre_alignmentCollection->getElementAt( j ) );

                                 // fixme has to be _cbcRealDUTsVec
                                 if (pre_alignment->getSensorID() == _cbcVirtualDUTId)
                                 {
                                        _dut_pre_align_x_1[i] = pre_alignment->getXOffset();
                                        _dut_pre_align_y_1[i] = pre_alignment->getYOffset();
                                        _dut_pre_align_z_1[i] = pre_alignment->getZOffset();
                                        _dut_pre_align_a_1[i] = pre_alignment->getAlpha();
                                        _dut_pre_align_b_1[i] = pre_alignment->getBeta();
                                        _dut_pre_align_c_1[i] = pre_alignment->getGamma();
                                 }

                                 // fixme has to be _cbcRealDUTsVec
                                 if (pre_alignment->getSensorID() == _cbcVirtualDUTId)
                                 {
                                        _dut_pre_align_x_2[i] = pre_alignment->getXOffset();
                                        _dut_pre_align_y_2[i] = pre_alignment->getYOffset();
                                        _dut_pre_align_z_2[i] = pre_alignment->getZOffset();
                                        _dut_pre_align_a_2[i] = pre_alignment->getAlpha();
                                        _dut_pre_align_b_2[i] = pre_alignment->getBeta();
                                        _dut_pre_align_c_2[i] = pre_alignment->getGamma();
                                 }
                        }
                } catch (lcio::DataNotAvailableException& e) {
                        streamlog_out( ERROR5 ) << "No pre_alignment collection with name " << _pre_alignmentCollectionName[i] << " !" << endl;
                }

                streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
                streamlog_out ( DEBUG7 ) << "Prealignment loading, iteration " << i+1 << ":" << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: X is:     " << _dut_pre_align_x_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: Y is:     " << _dut_pre_align_y_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: Z is:     " << _dut_pre_align_z_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: Alpha is: " << _dut_pre_align_a_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: Beta is:  " << _dut_pre_align_b_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT1 pre_alignment: Gamma is: " << _dut_pre_align_c_1[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: X is:     " << _dut_pre_align_x_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: Y is:     " << _dut_pre_align_y_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: Z is:     " << _dut_pre_align_z_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: Alpha is: " << _dut_pre_align_a_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: Beta is:  " << _dut_pre_align_b_2[i] << endl;
                streamlog_out ( DEBUG7 ) << "DUT2 pre_alignment: Gamma is: " << _dut_pre_align_c_2[i] << endl;
        }

        // get original prealignment
        if (_useOriginalPreAlignment == true)
        {

                _dut_original_pre_align_x_1 = 0.0;
                _dut_original_pre_align_y_1 = 0.0;
                _dut_original_pre_align_z_1 = 0.0;
                _dut_original_pre_align_a_1 = 0.0;
                _dut_original_pre_align_b_1 = 0.0;
                _dut_original_pre_align_c_1 = 0.0;

                _dut_original_pre_align_x_2 = 0.0;
                _dut_original_pre_align_y_2 = 0.0;
                _dut_original_pre_align_z_2 = 0.0;
                _dut_original_pre_align_a_2 = 0.0;
                _dut_original_pre_align_b_2 = 0.0;
                _dut_original_pre_align_c_2 = 0.0;

                try {
                        pre_alignmentCollection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _originalPreAlignment ) );

                        for (int j=0; j<_nTelPlanes;j++)
                        {
                                 EUTelAlignmentConstant * pre_alignment = static_cast< EUTelAlignmentConstant * > ( pre_alignmentCollection->getElementAt( j ) );

                                 // fixme has to be _cbcRealDUTsVec
                                 if (pre_alignment->getSensorID() == _cbcVirtualDUTId)
                                 {
                                        _dut_original_pre_align_x_1 = pre_alignment->getXOffset();
                                        _dut_original_pre_align_y_1 = pre_alignment->getYOffset();
                                        _dut_original_pre_align_z_1 = pre_alignment->getZOffset();
                                        _dut_original_pre_align_a_1 = pre_alignment->getAlpha();
                                        _dut_original_pre_align_b_1 = pre_alignment->getBeta();
                                        _dut_original_pre_align_c_1 = pre_alignment->getGamma();
                                 }
                                 // fixme has to be _cbcRealDUTsVec
                                 if (pre_alignment->getSensorID() == _cbcVirtualDUTId)
                                 {
                                        _dut_original_pre_align_x_2 = pre_alignment->getXOffset();
                                        _dut_original_pre_align_y_2 = pre_alignment->getYOffset();
                                        _dut_original_pre_align_z_2 = pre_alignment->getZOffset();
                                        _dut_original_pre_align_a_2 = pre_alignment->getAlpha();
                                        _dut_original_pre_align_b_2 = pre_alignment->getBeta();
                                        _dut_original_pre_align_c_2 = pre_alignment->getGamma();
                                 }
                        }
                } catch (lcio::DataNotAvailableException& e) {
                        streamlog_out( ERROR5 ) << "No original pre_alignment collection with name " << _originalPreAlignment << " !" << endl;
                }

                streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
                streamlog_out ( DEBUG7 ) << "Original Prealignment loading:" << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: X is:     " << _dut_original_pre_align_x_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: Y is:     " << _dut_original_pre_align_y_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: Z is:     " << _dut_original_pre_align_z_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: Alpha is: " << _dut_original_pre_align_a_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: Beta is:  " << _dut_original_pre_align_b_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_1 original pre_alignment: Gamma is: " << _dut_original_pre_align_c_1 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: X is:     " << _dut_original_pre_align_x_2 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: Y is:     " << _dut_original_pre_align_y_2 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: Z is:     " << _dut_original_pre_align_z_2 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: Alpha is: " << _dut_original_pre_align_a_2 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: Beta is:  " << _dut_original_pre_align_b_2 << endl;
                streamlog_out ( DEBUG7 ) << "DUT_2 original pre_alignment: Gamma is: " << _dut_original_pre_align_c_2 << endl;
        }
}

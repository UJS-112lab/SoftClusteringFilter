#if !defined(_SCF_H_INCLUDED_)
#define _SCF_H_INCLUDED_

using namespace std;

class SCF
{
public:
	SCF();
	virtual ~SCF();

	void PerformSCFO_ForGivenK(
		unsigned int*				ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		const int&					K,
		const double&				sigr);

private:
	void PerformSoftClusteringFiltering(
		unsigned int*				ubuff,
		double*&					kseedsr,
		double*&					kseedsg,
		double*&					kseedsb,
		double*&					kseedsx,
		double*&					kseedsy,
		const double&				sigr,
		const int&					STEP,
		const int&					numk,
		const int&					NUMITR);

	//============================================================================
	// Pick seeds when number of seeds is input.
	//============================================================================
	void GetRGBXYSeeds_ForGivenK(
		double*&					kseedsr,
		double*&					kseedsg,
		double*&					kseedsb,
		double*&					kseedsx,
		double*&					kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		double*&					edges,
		int*						seedIndices,
		int*						numseeds);

	//============================================================================
	// Move the seeds to low gradient positions to avoid putting seeds at region boundaries.
	//============================================================================
	void PerturbSeeds(
		double*&					kseedsr,
		double*&					kseedsg,
		double*&					kseedsb,
		double*&					kseedsx,
		double*&					kseedsy,
		double*&					edges,
		const int&					numseeds);

	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void DetectRgbEdges(
		const double*				rvec,
		const double*				gvec,
		const double*				bvec,
		const int&					width,
		const int&					height,
		double*&					edges);

private:
	int										m_width;
	int										m_height;
	int										m_depth;

	double*									m_rvec;
	double*									m_gvec;
	double*									m_bvec;

};

#endif // !defined(_SCF_H_INCLUDED_)

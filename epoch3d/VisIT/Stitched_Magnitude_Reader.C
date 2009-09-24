#include "Stitched_Magnitude_Reader.h"
#include "constants.h"
#include "CommonCommands.h"
#include <vtkDoubleArray.h>
#include <DebugStream.h>

#define min(a,b) (a<b?a:b)


//--------------------------------------------------------------------------------------------
 Stitched_Magnitude_Reader::Stitched_Magnitude_Reader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen,bool CacheOnly)
     : BlockReader(Handler,file,Owner,MaxStringLen,CacheOnly)
{
    //This is always a caching reader, so ignore CacheOnly
    //Note that despite this, it only caches the metadata, not the primary data
    
    //Can't guarantee that we'll be in the right place, so go there
    file->seekg(this->Owner->Offset,ios::beg);

    //Read the meshname/class
    this->MeshName=(char*)malloc(this->MaxStringLen+1);
    this->MeshClass=(char*)malloc(this->MaxStringLen+1);
    memset(this->MeshName,0,this->MaxStringLen+1);
    memset(this->MeshClass,0,this->MaxStringLen+1);
    file->read(this->MeshName,this->MaxStringLen);
    file->read(this->MeshClass,this->MaxStringLen);

    debug1 << "Stitched vector found on mesh " << MeshName<<endl;

    //Read dimensionality
    file->read((char*)&this->Dimensions,sizeof(int));
    SubBlocks=(Block**)malloc(this->Dimensions*sizeof(Block*));

}
//--------------------------------------------------------------------------------------------
Stitched_Magnitude_Reader::~Stitched_Magnitude_Reader()
{
    free(this->SubBlocks);
    free(this->MeshName);
    free(this->MeshClass);
}
//--------------------------------------------------------------------------------------------
void Stitched_Magnitude_Reader::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    char *Composite=(char*)malloc(2*MaxStringLen+1);
    memset(Composite,0,2*MaxStringLen+1);
    GetCompositeName(this->Owner->Name,this->Owner->Class,Composite);

    char *MeshComposite=(char*)malloc(2*MaxStringLen+1);
    memset(MeshComposite,0,2*MaxStringLen+1);
    GetCompositeName(this->MeshName,this->MeshClass,MeshComposite);

    debug1 << "Mesh associated with Block " << Owner->Name << " is " << MeshName << endl;

    avtCentering cent = AVT_NODECENT;

    avtScalarMetaData *smd = new avtScalarMetaData(Composite,MeshComposite,cent);
    md->Add(smd);

    free(Composite);
    free(MeshComposite);

}
//--------------------------------------------------------------------------------------------
vtkDataArray* Stitched_Magnitude_Reader::GetVar(int domain)
{
    //Spin past metadata
    file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length,ios::beg);

    char *Composite=(char*)malloc(2*MaxStringLen+1);
    char *Name=(char*)malloc(MaxStringLen+1);
    char *Class=(char*)malloc(MaxStringLen+1);

//Read the subblocks

    for (int i=0;i<this->Dimensions;++i)
    {
	memset(Composite,0,2*MaxStringLen+1);
	memset(Name,0,MaxStringLen+1);
	memset(Class,0,MaxStringLen+1);
	file->read(Name,this->MaxStringLen);
	file->read(Class,this->MaxStringLen);
	
	GetCompositeName(Name,Class,Composite);

	debug1 << "Searching for block " << Composite << endl;
	
	SubBlocks[i]=Handler->GetBlockByComposite(Composite);
	
	if (SubBlocks[i])
	{
	    if (SubBlocks[i]->Type != TYPE_MESH_VARIABLE)
	    {
		debug1 << "Subblock "<<i<<" is of incorrect type" <<endl;
		free(Composite);
		free(Name);
		free(Class);
		return NULL;
	    }
	    debug1 << "Subblock "<<i<<" is called "<<Name<<endl;
	}
	else
	{
	    debug1 << "Subblock "<<i << " missing" << endl;
	    free(Composite);
	    free(Name);
	    free(Class);
	    return NULL;
	}
    }

    free(Composite);
    free(Name);
    free(Class);

    //Don't support >3D data at the moment, will have to do some slicing when we do
    if (this->Dimensions > 3 || this->Dimensions <=0) return NULL;
    debug1 << "Ndims OK at " << this->Dimensions << endl;
    vtkDataArray *Data=NULL;

    bool DestroyReader;
    BlockReader *R=this->Handler->GetBlockReader(SubBlocks[0],DestroyReader);
    if (!R) return NULL;

    debug1 << "Asking for metadata from block " << SubBlocks[0]->Name << " of type " << SubBlocks[0]->Type << endl;

    Mesh_Variable_MetaData * MD=(Mesh_Variable_MetaData*)R->GetInternalMetaData();
    if (!MD || (this->Dimensions != MD->Dimensions && MD->Dimensions != DIMENSION_IRRELEVANT))
    {
	R->DestroyInternalMetaData(MD);
	if (DestroyReader) delete R;
	return NULL;
    }
    //There's much more data than this made available, but for the minute, ignore it
    this->SizeOfFloat=MD->SizeOfFloat;
    this->n_Elements=MD->n_Elements;

    //Finished with metadata
    R->DestroyInternalMetaData(MD);
    if (DestroyReader) delete R;

    if (this->SizeOfFloat==4) Data=vtkFloatArray::New();
    if (this->SizeOfFloat==8) Data=vtkDoubleArray::New();

    if (!Data) return NULL;

    Data->SetNumberOfTuples(this->n_Elements);

    long long nel_section=100051;
    void *v=malloc(this->SizeOfFloat*nel_section);
    void *pointdata=Data->GetVoidPointer(0);
    memset(pointdata,0,this->n_Elements*this->SizeOfFloat);
    void *pointnull=pointdata;
    for (int dimswing=0;dimswing <this->Dimensions;++dimswing)
    {
	file->seekg(SubBlocks[dimswing]->Offset + SubBlocks[dimswing]->Block_MD_Length,ios::beg);
	long long nel_left=this->n_Elements;
	while (nel_left >0) {
	    int imax=min(nel_section,nel_left);
	    file->read((char*)v,this->SizeOfFloat*imax);
	    if (this->SizeOfFloat == 4)
		for (int i=0;i<imax;++i)
		{
		    float *pointfloat=(float*)pointdata;
		    float readfloat=*(((float*)v)+i);
		    *pointfloat=(*pointfloat)+pow(readfloat,2.0);
                    pointdata=(void*)(((char*)pointdata)+this->SizeOfFloat);
		}
	    else
		for (int i=0;i<imax;++i)
		{
		    double *pointdouble=(double*)pointdata;
		    double readdouble=*(((double*)v)+i);
		    *pointdouble=(*pointdouble)+pow(readdouble,2.0);
                    pointdata=(void*)(((char*)pointdata)+this->SizeOfFloat);
		}
	    //Reduce remaining element count
	    nel_left=nel_left-nel_section;
	}
	pointdata=pointnull;
    }
    if (this->SizeOfFloat == 4)
	for (int iLoop=0;iLoop <this->n_Elements;++iLoop)
	{
	    float *pointfloat=(float*)pointdata;	
	    *pointfloat=pow(*pointfloat,0.5);
	    pointdata=(void*)(((char*)pointdata)+this->SizeOfFloat);
	}
    else
	for (int iLoop=0;iLoop <this->n_Elements;++iLoop)
	{
	    double *pointdouble=(double*)pointdata;	
	    *pointdouble=pow(*pointdouble,0.5);
	    pointdata=(void*)(((char*)pointdata)+this->SizeOfFloat);
	}


    free(v);
    
    return Data;
}
//--------------------------------------------------------------------------------------------

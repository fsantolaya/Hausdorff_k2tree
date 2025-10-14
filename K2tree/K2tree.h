/*
 * LKTree.h
 *  Es una extensión de los ktree que permite
 *  el etiquetado de aristas.
 *  la idea es poder etiquetar con una o más etiquetas una arista.
 *  cada valor de 32 bits.
 *  Created on: 15-11-2012
 *      Author: Miguel Romero Vásquez
 *      Editado: por Fernando Santolaya on 03/10/2025
 */
#ifndef LKTREE_H_
#define LKTREE_H_

#include <BitSequence.h>
using namespace cds_static;

#include "../Util/Point.h"

#ifndef uint
#define uint unsigned int
#endif

#ifndef K
#define K 2
#endif

typedef struct lnode {
	uint data;
	struct lnode** child;
} L_NODE;

typedef struct slkt {
	L_NODE * root;
	int max_Level;
	uint numberNodes;
	uint numberLeaves;
	uint numberTotalLeaves;
} lkt;

typedef struct matrixRep2{
	BitSequence * bt;		//Bitmap representando el arbol
	BitSequence * bn; 		//Bitmap representando el padre de las hojas
	BitSequence * bl;     	//Bitmap representando las hojas
	int maxLevel;			//Nivel maximo del arbol
	size_t numberOfNodes;
	size_t numberOfEdges;
	unsigned long long int  * div_level_table;
}MREP2;


L_NODE * createL_Node();
lkt * createLKTree(uint maxlevels);
int destroyLKTree(lkt * tree);
void insertNode(lkt * tree, int x, int y, uint label);
void _insertNode(lkt * tree, int x, int y, uint * labelArray);
int loadMREP2FromFile(MREP2 **ktree, FILE *file);
int saveMREP2ToFile(MREP2 *ktree, FILE *file);

MREP2 * createK2tree(string dataset1, int elevate, uint numberOfNodes);
void destroyK2tree(MREP2 * rep);
MREP2 * loadK2treeFromFile( const char *filename);
int saveK2treeToFile(MREP2 *snapshot, const char *filename) ;


/* AGREGADOS POR FERNANDO SANTOLAYA PARA EXPERIMENTOS */
double hausdorffDistHDK3MaxHeapv2(MREP2 * A, MREP2 * B);
double symmetricHausdorffDistance2D(MREP2 *A, MREP2 *B);
double naiveHDD( std::vector<Point>& p1,  std::vector<Point>& p2);
double hausdorffDistTaha2( vector<Point> &A,  vector<Point> &B);
extern int heapMaxElements;
extern int exitHeapElements;

int getHeapMaxElements();
int getHeapExitElements();

std::vector<Point> extractPointK2tree(MREP2 *J, int nelementos);
Point* extractPointK2tree2(MREP2 *J, int nelementos);

size_t sizeMREP2(MREP2 * rep);


#endif /* LKTREE_H_ */

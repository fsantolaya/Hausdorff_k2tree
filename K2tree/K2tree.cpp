/* *
 *  Created on: 15-11-2012
 *      Author: miguel
 *      Editado: Fernando Santolaya 03/10/2025
 */
#include "K2tree.h"
#include <libcdsBasics.h>
#include <limits>
#include <queue>
#include "../Util/Factory.h"

using namespace std;

int heapMaxElements = 0; // Definición
int exitHeapElements = 0; // Definición


//Operacioens con representación de árbol
//creación de la representación
L_NODE *createL_Node() {
    L_NODE *resp = (L_NODE *) malloc(sizeof(L_NODE));
    resp->data = 0;
    resp->child = NULL;
    return resp;
}
lkt *llenaK2tree(string ruta, int elevate) {
    lkt *tmp = createLKTree(elevate); //16 ERA ANTES
    string linea1;
    ifstream set1;
    set1.open(ruta);
    int z, i = 0;
    string token;
    int aux[2];
    while (getline(set1, linea1)) {
        //Ahora leemos la linea y obtenemos los puntos
        istringstream ss(linea1);
        z = 0;
        while (getline(ss, token, ',')) {
            aux[z] = atoi(token.c_str());
            z++;
        }
        insertNode(tmp, aux[0], aux[1], i);
        i++;
    }
    set1.close();

    return tmp;
}

lkt *createLKTree(uint maxlevels) {
    lkt *tree = (lkt *) malloc(sizeof(lkt));
    tree->root = createL_Node();
    tree->max_Level = maxlevels;
    tree->numberNodes = 0;
    tree->numberLeaves = 0;
    tree->numberTotalLeaves = 0;
    return tree;
}

//delete a sub-tree recursively
int destroyLNode(L_NODE *node) {
    int resp = 0;
    if (node == NULL) {
        return resp;
    } else {
        if (node->child != NULL) {
            //delete child recursively
            for (int i = 0; i < K * K; i++) {
                resp += destroyLNode(node->child[i]);
                node->child[i] = NULL;
            }
            free(node->child);
            node->child = NULL;
            resp += sizeof(L_NODE *) * K * K + sizeof(L_NODE *) * K * K;
        }
        free(node);
        node = NULL;
        resp += sizeof(L_NODE) + sizeof(L_NODE *);
    }
    return resp;
}

//elimina la reprecentación temporal de un lktree.
int destroyLKTree(lkt *tree) {
    int resp = 0;
    if (tree == NULL)
        return resp;
    //Elimino la raiz
    resp = destroyLNode(tree->root);
    tree->root = NULL;
    //elimino el nodo cabecera.
    free(tree);
    tree = NULL;
    return resp;
}

void insertNode(lkt *tree, int x, int y, uint label) {
    uint *pt = (uint *) malloc(sizeof(uint));
    *pt = label;
    _insertNode(tree, x, y, pt);
}

void _insertNode(lkt *tree, int x, int y, uint *labelArray) {
    uint i, node = 0;
    unsigned long long int div_level;
    int l = 0;
    L_NODE *n = tree->root;
    while (l <= tree->max_Level) {
        div_level = pow(K, tree->max_Level - l);
        node = (x / div_level) * K + y / div_level;
        if (l == tree->max_Level) {
            if (n->data == 0) {
                tree->numberLeaves++;
            }
            n->data = n->data | (0x1 << node); //1;
        } else {
            if (n->child == NULL) {
                if (l < tree->max_Level - 1) {
                    tree->numberNodes += K * K;
                } else {
                    tree->numberTotalLeaves += K * K;
                }
                n->child = (L_NODE **) calloc(K * K, sizeof(L_NODE *));
                for (i = 0; i < K * K; i++) {
                    n->child[i] = NULL;
                }
            }
            if (n->child[node] == NULL) {
                n->child[node] = createL_Node();
            }

            //n apunta al nodo insertado (intermedio)
            n = n->child[node];
        }
        x = x % div_level;
        y = y % div_level;
        l++;
    }
}

int saveMREP2ToFile(MREP2 *ktree, std::ofstream &ofs) {
    if (!ktree) return -1;

    // Guardar los campos del MREP2
    ofs.write(reinterpret_cast<const char *>(&ktree->maxLevel), sizeof(int));
    ofs.write(reinterpret_cast<const char *>(&ktree->numberOfNodes), sizeof(size_t));
    ofs.write(reinterpret_cast<const char *>(&ktree->numberOfEdges), sizeof(size_t));

    // Guardar el bitmap BT
    if (ktree->bt) {
        size_t btSize = ktree->bt->getSize();
        ofs.write(reinterpret_cast<const char *>(&btSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bt->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t btSize = 0;
        ofs.write(reinterpret_cast<const char *>(&btSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar el bitmap BN
    if (ktree->bn) {
        size_t bnSize = ktree->bn->getSize();
        ofs.write(reinterpret_cast<const char *>(&bnSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bn->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t bnSize = 0;
        ofs.write(reinterpret_cast<const char *>(&bnSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar el bitmap BL
    if (ktree->bl) {
        size_t blSize = ktree->bl->getSize();
        ofs.write(reinterpret_cast<const char *>(&blSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bl->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t blSize = 0;
        ofs.write(reinterpret_cast<const char *>(&blSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar la tabla de divisiones de nivel
    ofs.write(reinterpret_cast<const char *>(ktree->div_level_table),
              sizeof(unsigned long long int) * (ktree->maxLevel));

    return 0;
}

int loadMREP2FromFile(MREP2 **ktree, std::ifstream &ifs) {
    *ktree = (MREP2 *) malloc(sizeof(MREP2));
    if (!*ktree) {
        perror("Error al asignar memoria para MREP2");
        return -1;
    }

    // Leer los campos básicos
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->maxLevel), sizeof(int))) {
        perror("Error al leer maxLevel");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->numberOfNodes), sizeof(size_t))) {
        perror("Error al leer numberOfNodes");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->numberOfEdges), sizeof(size_t))) {
        perror("Error al leer numberOfEdges");
        free(*ktree);
        return -1;
    }

    // Cargar el bitmap BT
    size_t btSize;
    if (!ifs.read(reinterpret_cast<char *>(&btSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BT");
        free(*ktree);
        return -1;
    }
    if (btSize > 0) {
        (*ktree)->bt = BitSequence::load(ifs); // Asegúrate de que esta función maneje la memoria correctamente
    } else {
        (*ktree)->bt = nullptr;
    }

    // Cargar el bitmap BN
    size_t bnSize;
    if (!ifs.read(reinterpret_cast<char *>(&bnSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BN");
        free(*ktree);
        return -1;
    }
    if (bnSize > 0) {
        (*ktree)->bn = BitSequence::load(ifs);
    } else {
        (*ktree)->bn = nullptr;
    }

    // Cargar el bitmap BL
    size_t blSize;
    if (!ifs.read(reinterpret_cast<char *>(&blSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BL");
        free(*ktree);
        return -1;
    }
    if (blSize > 0) {
        (*ktree)->bl = BitSequence::load(ifs);
    } else {
        (*ktree)->bl = nullptr;
    }

    // Cargar la tabla de divisiones de nivel
    (*ktree)->div_level_table = (unsigned long long int *)
            malloc(sizeof(unsigned long long int) * ((*ktree)->maxLevel));
    if (!(*ktree)->div_level_table) {
        perror("Error al asignar memoria para div_level_table");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>((*ktree)->div_level_table),
                  sizeof(unsigned long long int) * ((*ktree)->maxLevel))) {
        perror("Error al leer div_level_table");
        free((*ktree)->div_level_table);
        free(*ktree);
        return -1;
    }

    return 0;
}


//Operaciones con representación compacta

MREP2 *createSnapshot(string dataset1, int elevate, uint cantPuntos) {
    //NO destruye lkt, hay que hacerlo desde afuera.
    //Estrategia general:
    //usando el arbol se construye el ktree de los puntos, y se construye un
    //arbol para las etiquetas.  Luego se construye el ktree para las estiquetas
    //con dicho arbol.

    //en este caso los valores de Numero de objetos y numero de arcos no se
    //si son necesarios, hay que revisar.
    int tamMatrix = pow(2, elevate);
    lkt *tree =  llenaK2tree(dataset1, elevate);
    MREP2 *ktree = (MREP2 *) malloc(sizeof(MREP2));
    ktree->maxLevel = tree->max_Level;
    ktree->numberOfNodes = tamMatrix;
    ktree->numberOfEdges = cantPuntos;
    size_t bits_BT_len = tree->numberNodes;
    size_t bits_BN_len = tree->numberTotalLeaves;
    size_t bits_LI_len = tree->numberLeaves * K * K;

    uint *bits_BT = new uint[uint_len(bits_BT_len, 1)]();
    uint *bits_BN = new uint[uint_len(bits_BN_len, 1)]();
    uint *bits_LI = new uint[uint_len(bits_LI_len, 1)]();

    uint k, j, queuecont, conttmp, node, pos = 0;

    unsigned long long div_level;
    int i;
    char isroot = 1;
    //Queue *q = createEmptyQueue();
    std::queue<L_NODE*> q;
    L_NODE *subTree;
    q.push(tree->root);
    queuecont = 1;

    for (i = 0; i < tree->max_Level; i++) {
        //i:iteración por nivel hasta el penúltimo
        conttmp = 0;
        div_level = pow(K, tree->max_Level - i);
        for (k = 0; k < queuecont; k++) {
            //por elemento de la cola de un mismo nivel
            subTree = (L_NODE *) q.front();
            q.pop();
            if ((subTree != NULL) && (subTree->child != NULL)) {
                for (j = 0; j < K * K; j++) {
                    node = j;
                    conttmp++;
                    q.push( subTree->child[node]);
                }
                if (!isroot)
                    bitset(bits_BT, pos);
            }
            if (!isroot)
                pos++;
            isroot = 0;
        }
        queuecont = conttmp;
    }

    ktree->bt = PlainBitSequenceFactory(bits_BT, bits_BT_len);
    pos = 0;
    size_t pos_inf = 0;

    //	fprintf(stderr,"Empezando bitmap de hojas utiles\n");
    while (!q.empty()) {
        subTree = (L_NODE *) q.front();
        q.pop();
        if ((subTree != NULL) && (subTree->data)) {
            //si no es 0
            bitset(bits_BN, pos);
            //hay que separar las listas

            for (i = 0; i < K * K; i++) {
                if ((subTree->data) & (0x1 << i)) {
                    bitset(bits_LI, pos_inf);
                }
                pos_inf++;
            }
        }
        pos++;
    }
   // destroyQueue(q);
    ktree->bn = PlainBitSequenceFactory(bits_BN, bits_BN_len);
    ktree->bl = PlainBitSequenceFactory(bits_LI, bits_LI_len);

    ktree->div_level_table = (long long uint *) malloc(
        sizeof(long long uint) * (ktree->maxLevel + 1));
    for (i = 0; i <= ktree->maxLevel; i++)
        ktree->div_level_table[i] = pow(K, ktree->maxLevel - i);

    delete[] bits_BT;
    bits_BT = NULL;
    delete[] bits_BN;
    bits_BN = NULL;
    delete[] bits_LI;
    bits_LI = NULL;
    return ktree;
}

void destroySnapshot(MREP2 *rep) {
    //destroyRepresentation(rep->ktree);
    delete rep->bl;
    delete rep->bn;
    delete rep->bt;
    free(rep->div_level_table);
    rep->bl = NULL;
    rep->bn = NULL;
    rep->bt = NULL;
    rep->div_level_table = NULL;
    free(rep);
    rep = NULL;
}

MREP2* loadSnapshotFromFile(const char *filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        perror("Error abriendo el archivo para cargar");
        return nullptr;
    }

    MREP2 * mrep = (MREP2 *) malloc(sizeof(MREP2));
    if (!mrep) {
        perror("Error al asignar memoria para MREP2");
        return nullptr;
    }

    // Leer el MREP2
    if (loadMREP2FromFile(&mrep, ifs) != 0) {
        perror("Error al cargar MREP2 desde archivo");
        free(mrep);
        return nullptr;
    }

    return mrep;
}

// Guardar Snapshot en archivo
int saveSnapshotToFile(MREP2 *snapshot, const char *filename) {
    if (!snapshot || !filename) return -1;

    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        perror("Error abriendo el archivo para guardar el Snapshot");
        return -1;
    }
    // Guardar el MREP2 (árbol k2tree)
    if (saveMREP2ToFile(snapshot, ofs) != 0) {
        ofs.close();
        return -1;
    }
    ofs.close();
    return 0;
}


void firstLeaf(MREP2 *rep, uint p1, uint p2, uint q1, uint q2,
               uint dp, uint dq, int x, int l, uint &n, Point &p) {
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //Point p;
                    p.setX(dp + i);
                    p.setY(dq + j);
                    return;
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                    //snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,dq + j, Oid, X, Y, n);
                }
            }
            leaf += K;
        }
    }

    if ((l == rep->maxLevel - 1)
        && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                firstLeaf(rep, 0, 0, 0, 0, dp + i, dq + j,
                          y + K * i + j, l + 1, n, p);
            }
        }
    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                firstLeaf(rep, p1new, p2new, q1new, q2new,
                          dp + divlevel * i, dq + divlevel * j, y + K * i + j,
                          l + 1, n, p);
            }
        }
    }
}

typedef struct sArea {
    unsigned long int x1, x2, y1, y2;
} SpatialArea;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * element Queue
 */
typedef struct sElementQueue {
    long int sTree;
    SpatialArea area;
    double priority; //la distancia...
} ElementQueue;

typedef struct sElementQueueSHD {
    long int sTree;
    SpatialArea area;
    char direccion;
    double priority; //la distancia...
} ElementQueueSHD;

void printElementQueue(ElementQueue e) {
    printf("ElementQueue:\n");
    printf("  sTree: %ld\n", e.sTree);
    printf("  Area: [x1=%lu, x2=%lu, y1=%lu, y2=%lu]\n",
           e.area.x1, e.area.x2, e.area.y1, e.area.y2);
    printf("  Priority (distancia): %.6f\n", e.priority);
}

typedef struct sCandidateQueue {
    long int oid;
    long int x;
    long int y;
    double priority;
} CandidateElement;


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Max Heap comparator
 */
struct maxHeapComparator {
    bool operator()(CandidateElement i, CandidateElement j) {
        return i.priority < j.priority;
    }
};

struct maxHeapComparatorHD {
    bool operator()(ElementQueue i, ElementQueue j) {
        return i.priority < j.priority;
    }
};

struct maxHeapComparatorSHD {
    bool operator()(ElementQueueSHD i, ElementQueueSHD j) {
        return i.priority < j.priority;
    }
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Min Heap comparator
 */
struct minHeapComparator {
    bool operator()(ElementQueue i, ElementQueue j) {
        return i.priority > j.priority;
    }
};

struct minHeapComparatorSHD {
    bool operator()(ElementQueueSHD i, ElementQueueSHD j) {
        return i.priority > j.priority;
    }
};

typedef struct sPriorityQueue {
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pPila;
    int bandera;
} PriorityQueue;

typedef struct sHdQueue {
    PriorityQueue pPila;
    ElementQueue element;
    int bandera;
} HdQueue;

double maxDist(Point p, SpatialArea a) {
    long int dx1 = (p.getX() - a.x1);
    dx1 = dx1 >= 0 ? dx1 : -dx1;

    long int dx2 = (p.getX() - a.x2);
    dx2 = dx2 >= 0 ? dx2 : -dx2;

    long int dy1 = p.getY() - a.y1;
    dy1 = dy1 >= 0 ? dy1 : -dy1;

    long int dy2 = p.getY() - a.y2;
    dy2 = dy2 >= 0 ? dy2 : -dy2;

    long int dx = max(dx1, dx2);
    long int dy = max(dy1, dy2);
    return sqrt((dx * dx) + (dy * dy));
    //return (dx * dx) + (dy * dy);
}
double minDist(Point p, SpatialArea a) {
    long int x = p.getX();
    long int y = p.getY();
    if (x < a.x1) {
        //hay 3 casos o bien está dentro del rango y1,y2 o fuera.
        if (y < a.y1) {
            long int dx = a.x1 - x;
            long int dy = a.y1 - y;
            return sqrt((dx * dx) + (dy * dy));
            //return (dx * dx) + (dy * dy);
        } else if (y <= a.y2) {
            return double(a.x1 - x);
        } else {
            long int dx = a.x1 - x;
            long int dy = y - a.y2;
            return sqrt((dx * dx) + (dy * dy));
            //return (dx * dx) + (dy * dy);
        }
    } else if (x <= a.x2) {
        if (y < a.y1) {
            return double(a.y1 - y);
        } else if (y <= a.y2) {
            return min(min((x - a.x1), (a.x2 - x)), min((y - a.y1), (a.y2 - y)));
        } else {
            return double(y - a.y2);
        }
    } else {
        if (y < a.y1) {
            long int dx = x - a.x2;
            long int dy = a.y1 - y;
            return sqrt((dx * dx) + (dy * dy));
            //return (dx * dx) + (dy * dy);
        } else if (y <= a.y2) {
            return double(x - a.x2);
        } else {
            long int dx = x - a.x2;
            long int dy = y - a.y2;
            return sqrt((dx * dx) + (dy * dy));
            //return (dx * dx) + (dy * dy);
        }
    }
    return 0.0;
}
double maxDistSinSqrt(Point p, SpatialArea a) {
    long int dx1 = (p.getX() - a.x1);
    dx1 = dx1 >= 0 ? dx1 : -dx1;

    long int dx2 = (p.getX() - a.x2);
    dx2 = dx2 >= 0 ? dx2 : -dx2;

    long int dy1 = p.getY() - a.y1;
    dy1 = dy1 >= 0 ? dy1 : -dy1;

    long int dy2 = p.getY() - a.y2;
    dy2 = dy2 >= 0 ? dy2 : -dy2;

    long int dx = max(dx1, dx2);
    long int dy = max(dy1, dy2);
    //return sqrt((dx * dx) + (dy * dy));

    return (dx * dx) + (dy * dy);
}
double minDistSinSqrt(Point p, SpatialArea a) {
    long int x = p.getX();
    long int y = p.getY();
    if (x < a.x1) {
        //hay 3 casos o bien está dentro del rango y1,y2 o fuera.
        if (y < a.y1) {
            long int dx = a.x1 - x;
            long int dy = a.y1 - y;

            return (dx * dx) + (dy * dy);
            //return (dx * dx) + (dy * dy);
        } else if (y <= a.y2) {
            return double((a.x1 - x) * (a.x1 - x));
        } else {
            long int dx = a.x1 - x;
            long int dy = y - a.y2;
            return (dx * dx) + (dy * dy);
            //return (dx * dx) + (dy * dy);
        }
    } else if (x <= a.x2) {
        if (y < a.y1) {
            return double((a.y1 - y) * (a.y1 - y));
        } else if (y <= a.y2) {
            //return min(min((x - a.x1), (a.x2 - x)), min((y - a.y1), (a.y2 - y)));
            double aux = min(min((x - a.x1), (a.x2 - x)), min((y - a.y1), (a.y2 - y)));
            return aux * aux;
        } else {
            return double((y - a.y2) * (y - a.y2));
        }
    } else {
        if (y < a.y1) {
            long int dx = x - a.x2;
            long int dy = a.y1 - y;

            return (dx * dx) + (dy * dy);
            //return (dx * dx) + (dy * dy);
        } else if (y <= a.y2) {
            return double((x - a.x2) * (x - a.x2));
        } else {
            long int dx = x - a.x2;
            long int dy = y - a.y2;

            //return sqrt((dx * dx) + (dy * dy));
            return (dx * dx) + (dy * dy);
        }
    }

    return 0.0;
}

//Funcion que calcula la mayor distancia entre areas
double maxDistAreas(SpatialArea s, SpatialArea b) {
    Point p;
    double d, d2, d3, d4;
    //Esto debemos optimizarlo
    //Obtenemos el primer vertice
    p.setX(s.x1);
    p.setY(s.y1);
    d = maxDistSinSqrt(p, b);
    //Obtenemos el segundo vertice
    p.setX(s.x1);
    p.setY(s.y2);
    d2 = maxDistSinSqrt(p, b);
    if (d2 > d) {
        d = d2;
    }
    //Obtenemos el tercer vertice
    p.setX(s.x2);
    p.setY(s.y1);
    d3 = maxDistSinSqrt(p, b);
    if (d3 > d) {
        d = d3;
    }
    //Obtenemos el cuarto vertice
    p.setX(s.x2);
    p.setY(s.y2);
    d4 = maxDistSinSqrt(p, b);
    if (d4 > d) {
        d = d4;
    }
    return d;
}

inline int isLeaf(MREP2 *mrep, long int index) {
    //index está mapeado como si la implementación fuera un único bitmap
    //con la finalidad de encapsular las diferentes implementaciones existentes
    //del k^2-tree
    //verdadero si es una hoja de último nivel sino falso.
    //si index es -1 entonces es la raíz
    return (index != -1)
           && ((mrep->bt->getLength() + mrep->bn->getLength()) <= index);
}
inline int posInLeaf(MREP2 *mrep, long int index) {
    return index - (mrep->bt->getLength() + mrep->bn->getLength());
}
inline bool esCero(MREP2 *mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    if (x == -1) {
        return false;
    } else if (x < btLen) {
        return mrep->bt->access(x) == 0;
    } else if (x < (btLen + bnLen)) {
        return mrep->bn->access(x - btLen) == 0;
    } else {
        //es una hoja de último nivel.
        return mrep->bl->access(x - (btLen + bnLen)) == 0;
    }
}

inline long int firstChild(MREP2 *mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    long int btOnes = mrep->bt->countOnes();
    if (x == -1) {
        return 0;
    } else if (x < btLen) {
        return (mrep->bt->access(x) == 0) ? -1 : mrep->bt->rank1(x) * K * K;
    } else if (x < btLen + bnLen) {
        /*		cout << endl << "btOnes:" << btOnes << " x in bt:" << x - btLen
         << " rank1 de x en bt:" << mrep->bn->rank1(x - btLen) << endl;
         cout << ((btOnes + (mrep->bn->rank1(x - btLen))) * K * K) << endl;
         */
        return (mrep->bn->access(x - btLen) == 0)
                   ?
                   //error: -1 : (btLen + (mrep->bn->rank1(x - btLen) * K * K));
                   -1
                   : ((btOnes + (mrep->bn->rank1(x - btLen))) * K * K);
    } else {
        //es una hoja de último nivel.
        return -1;
    }
}

//Esto está implementada para un K=2, pero hay que hacerlo para cualquier valor de K.
inline SpatialArea getSubArea(SpatialArea a, int i) {
    SpatialArea resp;
    switch (i) {
        case 0:
            resp.x1 = a.x1;
            resp.x2 = (a.x1 + a.x2) / 2;
            resp.y1 = a.y1;
            resp.y2 = (a.y1 + a.y2) / 2;
            break;
        case 1:
            resp.x1 = a.x1;
            resp.x2 = (a.x1 + a.x2) / 2;
            resp.y1 = ((a.y1 + a.y2) / 2) + 1;
            resp.y2 = a.y2;

            break;
        case 2:
            resp.x1 = ((a.x1 + a.x2) / 2) + 1;
            resp.x2 = a.x2;
            resp.y1 = a.y1;
            resp.y2 = (a.y1 + a.y2) / 2;

            break;
        case 3:
            resp.x1 = ((a.x1 + a.x2) / 2) + 1;
            resp.x2 = a.x2;
            resp.y1 = ((a.y1 + a.y2) / 2) + 1;
            resp.y2 = a.y2;
            break;
    }
    return resp;
}


/**************************************************************/
void nnMax(MREP2 *rep, Point q, double &max) {
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
    double c = std::numeric_limits<double>::max(); //
    ElementQueue e;
    e.sTree = -1;
    e.area.x1 = 0;
    //el 2 es porque k = 2, reemplazar por K.
    e.area.x2 = pow(K, rep->maxLevel + 1) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.priority = maxDistSinSqrt(q, e.area); //Se calcula la maxima distancia hacia los vertices
    long int subTreeIndex;
    double minDis2Child;
    pQueue.push(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        if (isLeaf(rep, e.sTree)) {
            //si es hoja de último nivel
            if ((e.priority <= max)) {
                //return max;
                return;
            } else {
                if (e.priority < c) {
                    c = e.priority;
                }
            }
            //cout << "HOJA?: " << e.priority << endl;
        } else {
            subTreeIndex = firstChild(rep, e.sTree);
            //child function return -1 if sTree no have any child
            //no tiene hijo si en la posición e.sTreee hay un 0.

            for (int i = 0; i < (K * K); i++) {
                if (!esCero(rep, subTreeIndex + i)) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);
                    child.priority = maxDistSinSqrt(q, child.area);
                    if (child.priority <= max) {
                        //return max;
                        return;
                    }
                    minDis2Child = minDistSinSqrt(q, child.area);
                    if (minDis2Child < c) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
    max = c;
    //cout << "CMAX: " << sqrt(max) << endl;
}


/* agregados por Fernando Santolaya para rastrear tamaños del heap */
void saveMaxPqueue(unsigned long num) {
    if (num > heapMaxElements) {
        heapMaxElements = num;
    }
}

void saveExitHeap(unsigned long num) {
    exitHeapElements = num;
}

////////////////////////////////////////
//-------------------------------------------------------------------------------------------------------------------------------
// D I S T A N C I A  S I M E T R I C A  D E  H A U S D O R F  B A S A D O  E N   H D K P 2
//------------------------------------------------------------------------------------------------------------------------------
// Calcula la distancia simétrica de hausdorff.
// este código está en estado inicial. No usar, porque no calcula la simetrica aún.
//////////////////////////////////////////////////////////////////
/// AQUI MODIFICO PARA CREAR LA VERSIÓN QUE CALCULA HDSIM.
//este es el hausdorff detonado con caleta de poda!!!!
//VERSION 2D

//Aqui se busca el punto más cercano en "B" de un área de "A". Si el sentido es de B->A,
//cambian los roles en la llamada, aunque internamente se llama B. retorna -1 si no se
//puede mejorar el supermax (early break).

double isCandidateSHD(const ElementQueueSHD &eq, MREP2 *B, double &supermax) {
    priority_queue<ElementQueueSHD, std::vector<ElementQueueSHD>, minHeapComparatorSHD> pQueue;
    long int subTreeIndex;
    ElementQueueSHD e;
    e.sTree = -1;
    e.area.x1 = 0;
    //el 2 es porque k = 2, reemplazar por K.
    e.area.x2 = pow(K, B->maxLevel + 1) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    pQueue.push(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        if (isLeaf(B, e.sTree)) {
            //return true;
            return e.priority;
        } else {
            subTreeIndex = firstChild(B, e.sTree);
            for (int i = 0; i < (K * K); i++) {
                if (!esCero(B, subTreeIndex + i)) {
                    ElementQueueSHD child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);
                    child.priority = maxDistAreas(eq.area, child.area);
                    if (child.priority <= supermax) {
                        //return false;
                        return -1;
                    }
                    pQueue.push(child);
                }
            }
        }
    }
    //significa que B está vacío.
    return -1;
}

//-------------------------invento de miguel romero..:) dont work..---------------------------
void symmetricHausdorffDistance2D(MREP2 *A, MREP2 *B, double &supermax) {
    priority_queue<ElementQueueSHD, std::vector<ElementQueueSHD>, maxHeapComparatorSHD> pQueue;
    //Calculamos la region de A
    ElementQueueSHD aa;
    aa.direccion = '>'; //A-->B
    //inicialmente ambos cuadrantes son iguales, solo cambia la dirección
    aa.sTree = -1;
    aa.area.x1 = 0;
    aa.area.x2 = pow(K, A->maxLevel + 1) - 1;
    aa.area.y1 = aa.area.x1;
    aa.area.y2 = aa.area.x2;
    aa.priority = A->numberOfNodes * sqrt(2); // la diagonal es la mayor distancia posible..
    pQueue.push(aa);

    //Calculamos la region de B
    ElementQueueSHD bb;
    bb.direccion = '<'; //B-->A
    //inicialmente ambos cuadrantes son iguales, solo cambia la dirección
    bb.sTree = -1;
    bb.area.x1 = 0;
    bb.area.x2 = pow(K, B->maxLevel + 1) - 1;
    bb.area.y1 = bb.area.x1;
    bb.area.y2 = bb.area.x2;
    bb.priority = B->numberOfNodes * sqrt(2); // la diagonal es la mayor distancia posible..
    pQueue.push(bb);

    long int subTreeIndex;

    while (!pQueue.empty()) {
        // para saber el maximo ocupado por pqueue y reportar el espacio ocupado..
        saveMaxPqueue(pQueue.size());
        aa = pQueue.top();
        //hay que reemplazar toda ocurrencia de A con k2treeA
        //que será A o B, según corresponda a la dirección
        MREP2 *k2treeA = aa.direccion == '>' ? A : B;
        MREP2 *k2treeB = aa.direccion == '>' ? B : A;
        // y la fhecla no va al reves? No, porque es el B, si va > es A->B
        //'>' caso normal(K2treeA=A, K2treeB=B), sino al reves.
        pQueue.pop();
        /* Si la mayor distancia posible del cuadrante aa es menor que supermax
         * ningun cuadrante del heap (pqueue) mejorará el resultado conocido y por tanto se acaba la busqueda.
          * */
        if (aa.priority <= supermax) {
            saveExitHeap(pQueue.size()); //agregado por fernando.
            return;
        }
        //Aqui verificamos si es hoja.
        if (isLeaf(k2treeA, aa.sTree)) {
            //cout << "Hoja" << endl;
            Point q(aa.area.x1, aa.area.y1);
            double aux = supermax;
            nnMax(k2treeB, q, supermax);
            if ((aux >= supermax)) {
                //supermax es pasada por referencia, se actualiza sola
                supermax = aux; //Ahora bien, si no mejoramos supermax, debemos volverla a su valor anterior
            }
        } else {
            //hay que reemplazar toda ocurrencia de A con k2treeA
            //que será A o B, según corresponda a la dirección
            MREP2 *k2treeA = aa.direccion == '>' ? A : B;
            subTreeIndex = firstChild(k2treeA, aa.sTree);
            for (int i = 0; i < (K * K); i++) {
                if (!esCero(k2treeA, subTreeIndex + i)) {
                    ElementQueueSHD child;
                    //el hijo mantiene la dirección del cálculo del padre.
                    child.direccion = aa.direccion;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(aa.area, i);
                    MREP2 *k2treeB = aa.direccion == '>' ? B : A;
                    //TODO: Revisar bien este paso...
                    //isCandidateSHD indica -1 si no es candidato, en caso de serlo,
                    //refina la distancia del cuadrante child al punto más cercano y real
                    //en el k2tree correspondiente al sentido de la búsqueda.
                    if ((child.priority = isCandidateSHD(child, k2treeB, supermax)) != -1) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
    saveExitHeap(pQueue.size()); //agregado por fernando.

    //cout << "N° de Puntos MaxHeap: " << contador << endl;
}

/***************************************/
//////////////////////////////////////////////////////////////////
double isCandidate(ElementQueue eq, MREP2 *B, double &supermax) {
priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
long int subTreeIndex;
ElementQueue e;
e.sTree = -1;
e.area.x1 = 0;
//el 2 es porque k = 2, reemplazar por K.
e.area.x2 = pow(K, B->maxLevel + 1) - 1;
e.area.y1 = e.area.x1;
e.area.y2 = e.area.x2;
pQueue.push(e);
while (!pQueue.empty()) {
    e = pQueue.top();
    pQueue.pop();
    if (isLeaf(B, e.sTree)) {
        //return true;
        //return -1;
        return e.priority;
    } else {
        subTreeIndex = firstChild(B, e.sTree);
        for (int i = 0; i < (K * K); i++) {
            if (!esCero(B, subTreeIndex + i)) {
                ElementQueue child;
                child.sTree = subTreeIndex + i;
                child.area = getSubArea(e.area, i);
                child.priority = maxDistAreas(eq.area, child.area);
                if (child.priority <= supermax) {
                    //return child.priority;
                    //return false;
                    return -1;
                }
                pQueue.push(child);
            }
        }
    }
}
return -1;
}

//este es el hausdorff detonado con caleta de poda!!!!
void hdkMaxHeapv2(MREP2 *A, MREP2 *B, double &supermax) {
    priority_queue<ElementQueue, std::vector<ElementQueue>, maxHeapComparatorHD> pQueue;
    //Calculamos la region de A
    ElementQueue aa;
    aa.sTree = -1;
    aa.area.x1 = 0;
    //el 2 es porque k = 2, reemplazar por K.
    aa.area.x2 = pow(K, A->maxLevel + 1) - 1;
    aa.area.y1 = aa.area.x1;
    aa.area.y2 = aa.area.x2;
    aa.priority = A->numberOfNodes * sqrt(2); // la diagonal es la mayor distancia posible..
    int contador = 0;
    long int subTreeIndex;
    PriorityQueue pfilaB;
    pQueue.push(aa);
    while (!pQueue.empty()) {
        saveMaxPqueue(pQueue.size()); // para saber el maximo ocupado por pqueue y reportar el espacio ocupado..
        aa = pQueue.top();
        pQueue.pop();
        /* Si la mayor distancia posible del cuadrante aa es menor que supermax
         * ningun cuadrante del heap (pqueue) mejorará el resultado conocido y por tanto se acaba la busqueda.
          * */
        if (aa.priority <= supermax) {
            saveExitHeap(pQueue.size()); //agregado por fernando.
            return;
        }
        //Aqui verificamos si es hoja.
        if (isLeaf(A, aa.sTree)) {
            //cout << "Hoja" << endl;
            Point q;
            q.setX(aa.area.x1);
            q.setY(aa.area.y1);
            double aux = supermax;
            contador++;
            nnMax(B, q, supermax);
            if ((aux >= supermax)) {
                //supermax es pasada por referencia, se actualiza sola
                supermax = aux; //Ahora bien, si no mejoramos supermax, debemos volverla a su valor anterior
            }
        } else {
            subTreeIndex = firstChild(A, aa.sTree);
            for (int i = 0; i < (K * K); i++) {
                if (!esCero(A, subTreeIndex + i)) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(aa.area, i);
                    if ((child.priority = isCandidate(child, B, supermax)) != -1) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
    saveExitHeap(pQueue.size()); //agregado por fernando.

    //cout << "N° de Puntos MaxHeap: " << contador << endl;
}


/*AGREGADO POR FERNANDO SANTOLAYA PARA HACER EXPERIMENTOS NUEVOS*/

//esta es el ultimo hausdorff con caleta de poda
double hausdorffDistHDK3MaxHeapv2(MREP2 *A, MREP2 *B) {
    double supermax = 0.0;
    hdkMaxHeapv2(A, B, supermax);
    return sqrt(supermax);
}
//symetricHausdorff
double symmetricHausdorffDistance2D(MREP2 *A, MREP2 *B) {
    double supermax = 0.0;
    symmetricHausdorffDistance2D(A, B, supermax);
    return sqrt(supermax);
}


//*******************************************************************************************

double naiveHDD(Point p1[], Point p2[], int largo) {
    double d = 0.0;
    double cmin = 0.0;
    double cmax = 0.0;
    //for (Point x : p1) {
    //Point x = new Point();
    //Point y = new Point();
    for (int i = 0; i < largo; i++) {
        cmin = std::numeric_limits<double>::max();
        //for (Point y : p2) {
        for (int j = 0; j < largo; j++) {
            d = (p1[i].getX() - p2[j].getX()) * (p1[i].getX() - p2[j].getX()) +
                (p1[i].getY() - p2[j].getY()) * (p1[i].getY() - p2[j].getY());
            //distancia = sqrt((x1 - y1)*(x1 - y1 ) + (x2 - y2)*(x2 - y2));
            //d = p1[i].restar(p2[j]);
            //cout << d << ",";
            if (d < cmin) cmin = d;
        }
        //cout << endl;
        if (cmin > cmax) cmax = cmin;
    }
    delete[] p1;
    delete[] p2;
    return sqrt(cmax);
}

//entrada.
// rep: un snapshot
// q: el punto de referencia
//salida:
//	oid: el identificador del punto más cercano
//  p: la posición del punto más cercanco.
void NN(MREP2 *rep, Point q, uint &oid, Point &p) {
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
    ElementQueue e, c;
    e.sTree = 0;
    e.area.x1 = 0;
    //el 2 es porque k = 2, reemplazar por K.
    e.area.x2 = pow(K, rep->maxLevel) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.priority = maxDist(q, e.area);
    c = e;
    c.priority = e.priority + 1; //explicar el + 1 con un comentario!!!!
    int foundNN = 0;
    long int subTreeIndex;
    double minDis2Child;
    pQueue.push(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        if (isLeaf(rep, e.sTree)) {
            //si es hoja de último nivel
            if ((e.priority < c.priority)) {
                c = e;
                foundNN = 1;
                //Al encontrar un candidato, puede haber otro mejor
                //en otra área, hay que seguir buscando!
            }
        } else {
            subTreeIndex = firstChild(rep, e.sTree);
            //child function return -1 if sTree no have any child
            //no tiene hijo si en la posición e.sTreee hay un 0.
            if (subTreeIndex != -1) {
                for (int i = 0; i < (K * K); i++) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);
                    child.priority = maxDist(q, child.area);
                    minDis2Child = minDist(q, child.area);
                    if (minDis2Child < c.priority) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

int getHeapMaxElements() {
return heapMaxElements;
}

int getHeapExitElements() {
    return exitHeapElements;
}

void recursiveGetRank2(MREP2 *snap, uint p1, uint p2, uint q1, uint q2,
                      uint dp, uint dq, int x, int l, uint *&Oid, uint *&X,
                      uint *&Y, uint &n) {
    MREP2 *rep = snap;
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                   // snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,
                     //                        dq + j, Oid, X, Y, n);
                    X[n] = dp + i;
                    Y[n] = dq + j;
                    n++;
                }
            }
            leaf += K;
        }

    }

    if ((l == rep->maxLevel - 1)
        && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                recursiveGetRank2(snap, 0, 0, 0, 0, dp + i, dq + j,
                                 y + K * i + j, l + 1, Oid, X, Y, n);
            }
        }

    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                recursiveGetRank2(snap, p1new, p2new, q1new, q2new,
                                 dp + divlevel * i, dq + divlevel * j, y + K * i + j,
                                 l + 1, Oid, X, Y, n);
            }
        }

    }
}
void rangeQueryHDOriginal2(MREP2 *snap, uint p1, uint p2, uint q1, uint q2,
                          uint *&Oid, uint *&X, uint *&Y, uint &n) {
    Oid = new uint[snap->numberOfEdges]();
    X = new uint[snap->numberOfEdges]();
    Y = new uint[snap->numberOfEdges]();
    n = 0;
    recursiveGetRank2(snap, p1, p2, q1, q2, 0, 0, -1, -1, Oid, X, Y, n);
}
std::vector<Point> extractPointK2tree(MREP2 *J, int nelementos) {
    uint *Oid1;
    uint *X1;
    uint *Y1;
    uint m1 = pow(K, J->maxLevel + 1) - 1;
    uint n1;
    // Llama a la función rangeQueryHDOriginal para obtener los puntos
    rangeQueryHDOriginal2(J, 0, m1, 0, m1, Oid1, X1, Y1, n1);
    std::vector<Point> Ar(n1);
    // Llena el vector con los puntos obtenidos
    for (uint i = 0; i < n1; i++) {
        Point p1;
        p1.setX(X1[i]);
        p1.setY(Y1[i]);
        Ar[i] = p1;
    }
    // Retorna el vector de puntos
    return Ar;
}

size_t sizeMREP2(MREP2 *rep) {
    size_t totalByte = 0;
    if (rep != NULL) {
        totalByte = sizeof(MREP2);
        totalByte += sizeof(long long uint) * (rep->maxLevel + 1);
        totalByte += rep->bl->getSize();
        totalByte += rep->bn->getSize();
        totalByte += rep->bt->getSize();
    }
    return totalByte;
}



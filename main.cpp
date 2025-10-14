//
// Created by fernando santolaya on 09-09-24.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <K2tree.h>
#include <BitSequence.h>
#include "ryu-kamata.h"
#include "Util/TimeMesure.h"
#include "Util/utils.h"
#include "Util/MemoryMeasure.h"

using namespace std;
unsigned int nextPowerOfTwo(unsigned int n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return n + 1;
}
unsigned int matrizNecesaria(const std::string &filename) {
    std::ifstream set1(filename);
    if (!set1.is_open()) {
        throw std::runtime_error("No se pudo abrir el archivo");
    }

    std::string linea;
    std::string token;
    unsigned long long maxCoord = 0;

    while (std::getline(set1, linea)) {
        if (linea.empty()) continue;

        std::istringstream ss(linea);
        unsigned long long coords[2] = {0,0};
        int z = 0;

        while (std::getline(ss, token, ',')) {
            coords[z] = std::stoull(token);
            z++;
            if (z >= 2) break; // sólo x,y
        }

        if (z == 2) {
            maxCoord = std::max(maxCoord, std::max(coords[0], coords[1]));
        }
    }

    set1.close();

    // sumamos 1 porque si hay un punto en la coordenada maxCoord
    // la matriz debe tener espacio hasta ese índice
    return nextPowerOfTwo(maxCoord + 1);
}
std::string extractFilename(const std::string &path) {
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) {
        return path; // Path doesn't contain a '/' separator
    } else {
        return path.substr(pos + 1); // Extract filename from position after '/'
    }
}

std::string extractCategoryName(const std::string &path) {
    // Encuentra la última ocurrencia de "sets"
    size_t startPos = path.find("sets");
    if (startPos != std::string::npos) {
        startPos += 0; // Mueve el puntero después de "sets"

        // Encuentra la siguiente barra '/' después de "sets"
        size_t endPos = path.find('/', startPos);
        if (endPos != std::string::npos) {
            return path.substr(startPos, endPos - startPos);
        } else {
            // Si no se encuentra '/', toma el resto de la cadena
            return path.substr(startPos);
        }
    }
    return ""; // Retorna una cadena vacía si "sets" no se encuentra
}

std::string extractDistribution(const std::string &path) {
    size_t pos = path.find("Conjunto");
    if (pos != std::string::npos) {
        // Encuentra la posición del último número después de "Conjunto"
        pos = path.find_first_of("0123456789", pos);
        if (pos != std::string::npos && pos + 1 < path.length()) {
            char distChar = path[pos + 1]; // Captura el carácter justo después del número

            // Retorna el nombre completo de la distribución según la letra
            switch (distChar) {
                case 'R':
                    return "RANDOM";
                case 'M':
                    return "MIXTA";
                case 'G':
                    return "GAUSS";
                case 'T':
                    return "REAL";
                default:
                    return "UNKNOWN"; // En caso de que no sea R, M o G
            }
        }
    }
    return "UNKNOWN"; // Retorna "UNKNOWN" si no se encuentra la distribución
}

void printHeader() {
    printf("%-19s|%-8s|%-12s|%-12s|"
           "%-15s|%-15s|%-15s|%-15s|"
           "%-12s|%-7s|%-15s|%-15s|%-15s|%-15s|%-15s|\n",
           "SET NAME", "CATEGORY", "ALGORITHM", "DISTRIBUTION",
           "BUILD TIME", "QUERY TIME", "QUERY WALL TIME", "EXTRACT TIME",
           "Haus Dist", "LAMBDA", "STORAGE", "PEAK REAL MEM", "PEAK VIRT MEM",
           "MXHEAP", "EXIT HEAP"); // Añadir nuevos encabezados
}

void printStats(const Stats &stats) {
    printf("%-19s|%-8s|%-12s|%-12s|%15.01f|%15.01f|%15.01f|%15.01f|%12.01f|%7d|%15.01f|%15lu|%15lu|%15d|%15d|\n",
           stats.dataset.c_str(), stats.datasetCategory.c_str(),
           stats.algoritmo.c_str(), stats.distribution.c_str(),
           stats.buildTime, stats.queryTime, stats.queryTimeWall, stats.extractTime,
           stats.hauss_distance, stats.lambda, stats.storage,
           stats.peakRealMem, stats.peakVirtMem, stats.mxheap, stats.exit_heap); // Imprimir nuevos parámetros
}

int cuentaPuntos(string ruta1) {
    ifstream set1;
    set1.open(ruta1);
    string linea1;
    int nelementos1 = 0;
    //Necesitamos saber cuantos elementos tiene el
    while (getline(set1, linea1)) {
        nelementos1++;
    }
    //cout << nelementos1 << endl;
    set1.close();
    return nelementos1;
}

void imprimeResultados(Stats &resultados) {
    char *homedir = getenv("HOME");
    strcat(homedir, "/RESULTADOS/resultHauss.txt");
    if (checkFileExist(homedir) == 1 && fileIsEmpty(homedir) == 1) {
        printHeader();
        fflush(stdout);
    }
    printStats(resultados);
    fflush(stdout);
}

void K2TMAXHEAPv2(MREP2 *snapsp, MREP2 *snapsp2, Stats &resultados) {
    Cronometer *crono = cCronometer();
    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = std::max(hausdorffDistHDK3MaxHeapv2(snapsp, snapsp2),
                                  hausdorffDistHDK3MaxHeapv2(snapsp2, snapsp));
    time = (double) (stop_clock(crono) * 1000000.0);
    resultados.algoritmo = "K2T-MXHEAPv2";
    resultados.hauss_distance = haussCompleto;
    resultados.queryTime = time;
    resultados.queryTimeWall = wallTime(crono) * 1000000.0;
    resultados.storage = sizeMREP2(snapsp) + sizeMREP2(snapsp2);
    resultados.hauss_distance = haussCompleto;
    resultados.lambda = 0;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    resultados.exit_heap = getHeapExitElements();
    resultados.mxheap = getHeapMaxElements();
    imprimeResultados(resultados);
}

void K2TMAXHEAPv3(MREP2 *snapsp, MREP2 *snapsp2, Stats &resultados) {
    Cronometer *crono = cCronometer();
    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = symmetricHausdorffDistance2D(snapsp, snapsp2);
    time = (double) (stop_clock(crono) * 1000000.0);
    resultados.algoritmo = "symHD";
    resultados.hauss_distance = haussCompleto;
    resultados.queryTime = time;
    resultados.queryTimeWall = wallTime(crono) * 1000000.0;
    resultados.storage = sizeMREP2(snapsp) + sizeMREP2(snapsp2);
    resultados.hauss_distance = haussCompleto;
    resultados.lambda = 0;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    resultados.exit_heap = getHeapExitElements();
    resultados.mxheap = getHeapMaxElements();
    imprimeResultados(resultados);
}

void KAMATA(MREP2 *snapsp, MREP2 *snapsp2, Stats &resultados, int cantPuntosSet1, int cantPuntosSet2) {
    Cronometer *crono = cCronometer();
    double extrTime = 0.0;
    start_clock(crono);
    auto dataset1 = extractPointK2tree(snapsp, cantPuntosSet1);
    auto dataset2 = extractPointK2tree(snapsp2, cantPuntosSet2);
    extrTime = (double) (stop_clock(crono) * 1000000.0);

    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = hausKamata(dataset1, dataset2, 3);
    time = (double) (stop_clock(crono) * 1000000.0);

    resultados.algoritmo = "K2T-KAMATA";
    resultados.queryTime = time;
    resultados.queryTimeWall = wallTime(crono) * 1000000.0;
    resultados.extractTime = extrTime;
    resultados.storage = (((sizeof(int) * 2) * cantPuntosSet1) + ((sizeof(int) * 2) * cantPuntosSet2)) + (
                             sizeof(double) * cantPuntosSet1 + sizeof(double) * cantPuntosSet2);
    resultados.hauss_distance = haussCompleto;
    resultados.lambda = 3;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    imprimeResultados(resultados);
}

void NAIVE(MREP2 *snapsp, MREP2 *snapsp2, Stats &resultados, int cantPuntosSet1, int cantPuntosSet2) {
    Cronometer *crono = cCronometer();
    double extrTime = 0.0;
    start_clock(crono);
    auto dataset1 = extractPointK2tree(snapsp, cantPuntosSet1);
    auto dataset2 = extractPointK2tree(snapsp2, cantPuntosSet2);
    extrTime = (double) (stop_clock(crono) * 1000000.0);

    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = naiveHDD(dataset1, dataset2);
    time = (double) (stop_clock(crono) * 1000000.0);

    resultados.algoritmo = "K2T-NAIVE";
    resultados.queryTime = time;
    resultados.queryTimeWall = wallTime(crono) * 1000000.0;
    resultados.extractTime = extrTime;
    resultados.storage = (((sizeof(int) * 2) * cantPuntosSet1) + ((sizeof(int) * 2) * cantPuntosSet2)) + (
                             sizeof(double) * cantPuntosSet1 + sizeof(double) * cantPuntosSet2);
    resultados.hauss_distance = haussCompleto;
    resultados.lambda = 3;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    imprimeResultados(resultados);
}

void ejecutaExperimentos(const std::string &dataset1, const std::string &dataset2, int experimento) {
    Cronometer *crono = cCronometer();
    Stats resultados; // para guardar el resultado de los experimentos.

    int cantPuntosSet1 = cuentaPuntos(dataset1);
    int cantPuntosSet2 = cuentaPuntos(dataset2);

    double buildTime = 0.0;
    // int tamMatrix = pow(2, 16);

    // Aquí leerías el Snapshot desde el archivo
    std::string finalFilename1 = extractFilename(dataset1) + ".kt";
    std::string finalFilename2 = extractFilename(dataset2) + ".kt";

    MREP2 *snapsp = loadK2treeFromFile(finalFilename1.c_str());
    MREP2 *snapsp2 = loadK2treeFromFile(finalFilename2.c_str());

    resultados.dataset = extractFilename(dataset1);
    resultados.datasetCategory = extractCategoryName(dataset1);
    resultados.distribution = extractDistribution(resultados.dataset);
    resultados.buildTime = buildTime;

    // Ejecución del experimento según el tipo especificado
    switch (experimento) {
        case 1:
            KAMATA(snapsp, snapsp2, resultados, cantPuntosSet1, cantPuntosSet2);
            break;
        case 2:
            K2TMAXHEAPv2(snapsp, snapsp2, resultados);
            break;
        case 3: //EL NUEVO
            K2TMAXHEAPv3(snapsp, snapsp2, resultados);
            break;
        default:
            std::cerr << "Experimento desconocido." << std::endl;
            return;
    }

    // Liberar recursos
    destroyK2tree(snapsp); // Asegúrate de implementar esta función
    destroyK2tree(snapsp2); // Asegúrate de implementar esta función
    resultados.reset();
}

void test (const std::string &dataset1, const std::string &dataset2) {

    int cantPuntosSet1 = cuentaPuntos(dataset1);
    int cantPuntosSet2 = cuentaPuntos(dataset2);

    // Aquí leerías el Snapshot desde el archivo
    std::string finalFilename1 = extractFilename(dataset1) + ".kt";
    std::string finalFilename2 = extractFilename(dataset2) + ".kt";

    MREP2 *snapsp = loadK2treeFromFile(finalFilename1.c_str());
    MREP2 *snapsp2 = loadK2treeFromFile(finalFilename2.c_str());

    auto dataset11 = extractPointK2tree(snapsp, cantPuntosSet1);
    auto dataset22 = extractPointK2tree(snapsp2, cantPuntosSet2);

    cout <<"cantidad de puntos del set1: "<< dataset11.size() << endl;
    cout <<"cantidad de puntos del set2: "<< dataset22.size() << endl;

    auto k2treev3 = symmetricHausdorffDistance2D(snapsp, snapsp2);
    auto kamata = hausKamata(dataset11, dataset22, 3);
    auto taha = std::max(hausdorffDistTaha2(dataset11, dataset22), hausdorffDistTaha2(dataset22, dataset11));
    auto k2treev2 = std::max(hausdorffDistHDK3MaxHeapv2(snapsp, snapsp2),hausdorffDistHDK3MaxHeapv2(snapsp2, snapsp));
    auto naive =  std::max(naiveHDD(dataset11,dataset22),naiveHDD(dataset22,dataset11));
    std::cout <<"Hausdorff distance KAMATA = " << kamata << std::endl;
    std::cout <<"Hausdorff distance TAHA= " << taha << std::endl;
    std::cout <<"Hausdorff distance K2TREE V2= " << k2treev2 << std::endl;
    std::cout <<"Hausdorff distance K2TREE V3= " << k2treev3 << std::endl;
    std::cout <<"Hausdorff distance NAIVE= " << naive << std::endl;
    // Liberar recursos
    destroyK2tree(snapsp); // Asegúrate de implementar esta función
    destroyK2tree(snapsp2); // Asegúrate de implementar esta función
}

void construirEstructura(const std::string &dataset1, const std::string &dataset2) {
    MREP2 *snapsp = nullptr;
    MREP2 *snapsp2 = nullptr;

    unsigned int tamMatrix = matrizNecesaria(dataset1);
    std::cerr << "tam matrix is: " << tamMatrix << std::endl;
    unsigned int exponent = (unsigned int)std::log2(tamMatrix);
    cout <<"exponente is: "<< exponent<< endl;

    int cantPuntosSet1 = cuentaPuntos(dataset1);
    int cantPuntosSet2 = cuentaPuntos(dataset2);
    std::cerr << "cant point is is: " << cantPuntosSet1 << std::endl;
    std::cerr << "cant point2 is is: " << cantPuntosSet2 << std::endl;

    std::cerr << "creating k2trees" << std::endl;
    snapsp = createK2tree(dataset1,exponent, cantPuntosSet1);
    snapsp2 = createK2tree(dataset2, exponent,  cantPuntosSet2);

    std::string snap1filename = extractFilename(dataset1) + ".kt";
    std::string snap2filename = extractFilename(dataset2) + ".kt";
    std::cerr << "k2tree 1 filename is: " << snap1filename << std::endl;
    std::cerr << "k2tree 2 filename is: " << snap2filename << std::endl;

    std::cerr << "saving k2trees" << std::endl;
    saveK2treeToFile(snapsp, snap1filename.c_str());
    saveK2treeToFile(snapsp2, snap2filename.c_str());
}

int main(int argc, char **argv) {
    if (argc < 4) {
        printf("%s <PATH DATASETS> <EXPERIMENTO> <CONSTRUIR (1/0)>\n", argv[0]);
        printf("EXPERIMENTO: 1 (KAMATA), 2 (TAHA), 3 (K2T MX HEAP V2), 4 (K2T MX HEAP V3)\n");
        return -1;
    }

    std::string dataset1 = argv[1];
    std::string dataset2 = argv[2];

    // Debug prints
    std::cerr << "📂 Dataset 1: " << dataset1 << std::endl;
    std::cerr << "📂 Dataset 2: " << dataset2 << std::endl;


    int experimentoSeleccionado = std::stoi(argv[3]);
    int construir = std::stoi(argv[4]); // 1 para construir, 0 para no construir
    switch (construir) {
        case 1:
            construirEstructura(dataset1, dataset2);
            break;
        case 2:
            ejecutaExperimentos(dataset1, dataset2, experimentoSeleccionado);
            break;
        case 3:
            test(dataset1, dataset2);
            break;
    }

    return 0;
}

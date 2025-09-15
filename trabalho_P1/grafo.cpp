#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <queue>

class Grafo {
private:
    int numVertices;
    std::vector<std::vector<int>> adjList;

public:
    Grafo(int n) {
        numVertices = n;
        adjList.resize(n + 1);
    }

    void adicionarAresta(int u, int v) {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }

    static Grafo lerDeArquivo(const std::string& nomeArquivo) {
        std::ifstream arquivo(nomeArquivo);
        if (!arquivo.is_open()) {
            throw std::runtime_error("Erro ao abrir o arquivo!");
        }

        int n;
        arquivo >> n;
        Grafo g(n);

        int u, v;
        while (arquivo >> u >> v) {
            g.adicionarAresta(u, v);
        }

        arquivo.close();
        return g;
    }

    // ---------- Parte nova: Estatísticas ----------

    int contarArestas() const {
        int total = 0;
        for (int i = 1; i <= numVertices; i++) {
            total += adjList[i].size();
        }
        return total / 2; // divide por 2 porque o grafo é não-direcionado
    }

    std::vector<int> graus() const {
        std::vector<int> g;
        for (int i = 1; i <= numVertices; i++) {
            g.push_back(adjList[i].size());
        }
        return g;
    }

    // ---------- Componentes conexas ----------
    std::vector<std::vector<int>> componentesConexas() const {
        std::vector<bool> visitado(numVertices + 1, false);
        std::vector<std::vector<int>> componentes;

        for (int v = 1; v <= numVertices; v++) {
            if (!visitado[v]) {
                std::vector<int> componente;
                std::queue<int> fila;
                fila.push(v);
                visitado[v] = true;

                while (!fila.empty()) {
                    int u = fila.front();
                    fila.pop();
                    componente.push_back(u);

                    for (int viz : adjList[u]) {
                        if (!visitado[viz]) {
                            visitado[viz] = true;
                            fila.push(viz);
                        }
                    }
                }
                componentes.push_back(componente);
            }
        }
        return componentes;
    }

    // ---------- Gerar saída ----------
    void salvarEstatisticas(const std::string& nomeArquivo) {
        std::ofstream out(nomeArquivo);
        if (!out.is_open()) {
            throw std::runtime_error("Erro ao criar arquivo de saída!");
        }

        int nArestas = contarArestas();
        auto g = graus();

        int grauMin = *std::min_element(g.begin(), g.end());
        int grauMax = *std::max_element(g.begin(), g.end());
        double grauMedio = std::accumulate(g.begin(), g.end(), 0.0) / g.size();

        std::sort(g.begin(), g.end());
        double mediana;
        if (g.size() % 2 == 0) {
            mediana = (g[g.size()/2 - 1] + g[g.size()/2]) / 2.0;
        } else {
            mediana = g[g.size()/2];
        }

        out << "Numero de vertices: " << numVertices << "\n";
        out << "Numero de arestas: " << nArestas << "\n";
        out << "Grau minimo: " << grauMin << "\n";
        out << "Grau maximo: " << grauMax << "\n";
        out << "Grau medio: " << grauMedio << "\n";
        out << "Mediana do grau: " << mediana << "\n\n";

        // Componentes conexas
        auto comps = componentesConexas();
        out << "Numero de componentes conexas: " << comps.size() << "\n";
        for (size_t i = 0; i < comps.size(); i++) {
            out << "Componente " << i+1 << ": ";
            for (int v : comps[i]) {
                out << v << " ";
            }
            out << "\n";
        }

        out.close();
    }
};

// ----------------- Programa de teste -----------------
int main() {
    try {
        Grafo g = Grafo::lerDeArquivo("teste.txt");
        g.salvarEstatisticas("saida.txt");
        std::cout << "Arquivo de saida gerado com sucesso!\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }

    return 0;
}
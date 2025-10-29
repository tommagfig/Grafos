// grafo_refatorado.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <queue>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <string>
#include <stack>
#include <random>
#include <chrono>
#include <memory>
#include <cstdint>
#include <map>
using namespace std;

enum Representacao { LISTA, MATRIZ };

struct IGrafoRep {
    virtual ~IGrafoRep() = default;
    virtual void add_edge(int u, int v, double w) = 0;
    virtual void neighbors(int u, const function<void(int, double)>& f) const = 0;
    virtual int degree(int u) const = 0;
    virtual bool edge_exists(int u, int v) const = 0;
    virtual int n_vertices() const = 0;
};

// Nó da lista de adjacência
struct Node {
    int v;
    double w;
    Node* next;
    Node(int vert, double peso, Node* nxt = nullptr) : v(vert), w(peso), next(nxt) {}
};

class MinHeap {
private:
    vector<pair<int,double>> heap; // (vértice, distância)
    vector<int> pos;               // posição de cada vértice no heap (-1 se não estiver)

    void swapNodes(int i, int j) {
        swap(heap[i], heap[j]);
        pos[heap[i].first] = i;
        pos[heap[j].first] = j;
    }

    void siftUp(int i) {
        while (i > 0) {
            int pai = (i - 1) / 2;
            if (heap[pai].second <= heap[i].second) break;
            swapNodes(pai, i);
            i = pai;
        }
    }

    void siftDown(int i) {
        int size = heap.size();
        while (true) {
            int esq = 2*i + 1, dir = 2*i + 2, menor = i;
            if (esq < size && heap[esq].second < heap[menor].second) menor = esq;
            if (dir < size && heap[dir].second < heap[menor].second) menor = dir;
            if (menor == i) break;
            swapNodes(i, menor);
            i = menor;
        }
    }

public:
    explicit MinHeap(int n) : pos(n + 1, -1) {}

    bool empty() const { return heap.empty(); }

    void push(int v, double dist) {
        if (pos[v] != -1) { // já está no heap
            if (dist < heap[pos[v]].second) {
                heap[pos[v]].second = dist;
                siftUp(pos[v]);
            }
            return;
        }
        heap.push_back({v, dist});
        pos[v] = heap.size() - 1;
        siftUp(pos[v]);
    }

    pair<int,double> extract_min() {
        if (heap.empty()) return {-1, 1e18};
        auto minNode = heap[0];
        pos[minNode.first] = -1;
        if (heap.size() > 1) {
            heap[0] = heap.back();
            pos[heap[0].first] = 0;
        }
        heap.pop_back();
        if (!heap.empty()) siftDown(0);
        return minNode;
    }
};


//lista de adjacência
class ListaRep : public IGrafoRep {
private:
    vector<Node*> adj; 

public:
    explicit ListaRep(int n) : adj(n + 1, nullptr) {}

    void add_edge(int u, int v, double w) override {
        adj[u] = new Node(v, w, adj[u]);
        adj[v] = new Node(u, w, adj[v]);
    }

    void neighbors(int u, const function<void(int,double)>& f) const override {
        for (Node* curr = adj[u]; curr != nullptr; curr = curr->next) {
            f(curr->v, curr->w);
        }
    }

    int degree(int u) const override {
        int deg = 0;
        for (Node* curr = adj[u]; curr != nullptr; curr = curr->next) deg++;
        return deg;
    }

    bool edge_exists(int u, int v) const override {
        for (Node* curr = adj[u]; curr != nullptr; curr = curr->next)
            if (curr->v == v) return true;
        return false;
    }

    int n_vertices() const override {
        return (int)adj.size() - 1;
    }

    ~ListaRep() {
        for (auto head : adj) {
            Node* curr = head;
            while (curr) {
                Node* tmp = curr;
                curr = curr->next;
                delete tmp;
            }
        }
    }
};

//matriz de adjacência
class MatrizRep : public IGrafoRep {
private:
    vector<vector<double>> mat;
public:
    explicit MatrizRep(int n) : mat(n + 1, vector<double>(n + 1, 0.0)) {}
    void add_edge(int u, int v, double w) override {
        mat[u][v] = w;
        mat[v][u] = w;
    }
    void neighbors(int u, const function<void(int, double)>& f) const override {
        int n = (int)mat[u].size() - 1;
        for (int v = 1; v <= n; ++v) {
            if (mat[u][v] != 0.0) f(v, mat[u][v]);
        }
    }
    int degree(int u) const override {
        int cnt = 0;
        int n = (int)mat[u].size() - 1;
        for (int v = 1; v <= n; ++v) if (mat[u][v] != 0.0) ++cnt;
        return cnt;
    }
    bool edge_exists(int u, int v) const override {
        return mat[u][v] != 0.0;
    }
    int n_vertices() const override { return (int)mat.size() - 1; }
};

// ---------- Classe Grafo que usa a representação por composição ----------
class Grafo {
private:
    int n; // vértices
    int m; // arestas 
    unique_ptr<IGrafoRep> rep; // polimórfico

public:
    Grafo(int vertices, Representacao r) : n(vertices), m(0) {
        if (r == LISTA) rep = make_unique<ListaRep>(vertices);
        else rep = make_unique<MatrizRep>(vertices);
    }

    void adicionarAresta(int u, int v, double w) {
        if (u < 1 || u > n || v < 1 || v > n) return;
        if (!rep->edge_exists(u,v)) {
            rep->add_edge(u, v, w);
            ++m;
        }
    }

    static Grafo lerDeArquivo(const string& nomeArquivo, Representacao r) {
        ifstream entrada(nomeArquivo);
        if (!entrada.is_open()) {
            cerr << "Erro ao abrir o arquivo." << endl;
            exit(1);
        }
        int vertices, u, v;
        double w;
        entrada >> vertices;
        Grafo g(vertices, r);
        while (entrada >> u >> v >> w) {
            g.adicionarAresta(u, v, w);
        }
        entrada.close();
        return g;
    }

    int contarArestas() const {
        long long soma = 0;
        for (int i = 1; i <= n; ++i) soma += rep->degree(i);
        return (int)(soma / 2);
    }

    vector<int> graus() const {
        vector<int> g(n + 1, 0);
        for (int i = 1; i <= n; ++i) g[i] = rep->degree(i);
        return g;
    }

    int distancia(int origem, int destino) const {
        if (origem < 1 || origem > n || destino < 1 || destino > n) return -1;
        vector<int> dist(n + 1, -1);
        queue<int> q;
        dist[origem] = 0;
        q.push(origem);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            rep->neighbors(u, [&](int v, double peso) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    q.push(v);
                }
            });
            if (dist[destino] != -1) return dist[destino];
        }
        return dist[destino];
    }

    int diametro() const {
        int diam = 0;
        for (int i = 1; i <= n; ++i) {
            vector<int> dist(n + 1, -1);
            queue<int> q;
            dist[i] = 0;
            q.push(i);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                rep->neighbors(u, [&](int v, double peso) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        diam = max(diam, dist[v]);
                        q.push(v);
                    }
                });
            }
        }
        return diam;
    }

    vector<vector<int>> componentesConexas() const {
        vector<char> visitado(n + 1, 0);
        vector<vector<int>> componentes;
        for (int i = 1; i <= n; ++i) {
            if (!visitado[i]) {
                vector<int> comp;
                queue<int> q;
                q.push(i);
                visitado[i] = 1;
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    comp.push_back(u);
                    rep->neighbors(u, [&](int v, double peso) {
                        if (!visitado[v]) { visitado[v] = 1; q.push(v); }
                    });
                }
                componentes.push_back(std::move(comp));
            }
        }
        sort(componentes.begin(), componentes.end(),
             [](const vector<int>& a, const vector<int>& b){ return a.size() > b.size(); });
        return componentes;
    }

    void salvarEstatisticas(const string& nomeArquivo) const {
        ofstream out(nomeArquivo);
        if (!out.is_open()) throw runtime_error("Erro ao criar arquivo de saída!");

        int nArestas = contarArestas();
        vector<int> gAll = graus();
        vector<int> g(gAll.begin() + 1, gAll.end());

        int grauMin = *min_element(g.begin(), g.end());
        int grauMax = *max_element(g.begin(), g.end());
        double grauMedio = accumulate(g.begin(), g.end(), 0.0) / (double)n;

        vector<int> ordenado = g;
        sort(ordenado.begin(), ordenado.end());
        double mediana;
        if (n % 2 == 0) mediana = (ordenado[n/2 - 1] + ordenado[n/2]) / 2.0;
        else mediana = ordenado[n/2];

        int diam = diametro();
        auto comps = componentesConexas();

        out << fixed << setprecision(2);
        out << "Numero de vertices: " << n << "\n";
        out << "Numero de arestas: " << nArestas << "\n";
        out << "Grau minimo: " << grauMin << "\n";
        out << "Grau maximo: " << grauMax << "\n";
        out << "Grau medio: " << grauMedio << "\n";
        out << "Mediana do grau: " << mediana << "\n";
        out << "Diametro do grafo: " << diam << "\n\n";

        out << "Numero de componentes conexas: " << comps.size() << "\n\n";
        for (size_t i = 0; i < comps.size(); ++i) {
            out << "Componente " << i+1 << " (tamanho " << comps[i].size() << "): ";
            for (int v : comps[i]) out << v << " ";
            out << "\n";
        }
        out.close();
    }

    // BFS que escreve pai e nivel em arquivo 
    void BFS(int origem, const string& arquivo) const {
        vector<int> nivel(n + 1, -1);
        vector<int> pai(n + 1, 0);
        queue<int> q;
        nivel[origem] = 0;
        pai[origem] = 0;
        q.push(origem);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            vector<int> vizinhos;
            rep->neighbors(u, [&](int v, double peso) { vizinhos.push_back(v); });
            sort(vizinhos.begin(), vizinhos.end());

            for (int v : vizinhos) {
                if (nivel[v] == -1) {
                    pai[v] = u;
                    nivel[v] = nivel[u] + 1;
                    q.push(v);
                }
            }
        }

        ofstream out(arquivo);
        out << "BFS - Arvore de busca a partir do vertice " << origem << "\n";
        out << "Vertice\tPai\tNivel\n";
        for (int i = 1; i <= n; ++i)
            out << i << "\t\t" << pai[i] << "\t" << nivel[i] << "\n";
        out.close();
    }


    // DFS iterativa que escreve pai e nivel em arquivo
    void DFS(int origem, const string& arquivo) const {
        vector<int> nivel(n + 1, -1);
        vector<int> pai(n + 1, 0);

        stack<tuple<int,int,int>> s;
        s.push({origem, 0, 0});

        while (!s.empty()) {
            auto t = s.top();
            s.pop();
            int u = std::get<0>(t);
            int p = std::get<1>(t);
            int d = std::get<2>(t);

            if (nivel[u] != -1) continue;

            pai[u] = p;
            nivel[u] = d;

            vector<int> vizinhos;
            rep->neighbors(u, [&](int v, double peso){ vizinhos.push_back(v); });

            sort(vizinhos.begin(), vizinhos.end());
            for (auto it = vizinhos.rbegin(); it != vizinhos.rend(); ++it) {
                int v = *it;
                if (nivel[v] == -1) {
                    s.push({v, u, d + 1});
                }
            }
        }

        ofstream out(arquivo);
        out << "DFS - Arvore de busca a partir do vertice " << origem << "\n";
        out << "Vertice\tPai\tNivel\n";
        for (int i = 1; i <= n; ++i) out << i << "\t\t" << pai[i] << "\t" << nivel[i] << "\n";
        out.close();
    }


    // diametro aproximado por duas BFS (heurística)
    int diametroAproximado() const {
        auto bfsMaisDistante = [&](int origem){
            vector<int> dist(n + 1, -1);
            queue<int> q;
            dist[origem] = 0;
            q.push(origem);
            int ultimo = origem;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                ultimo = u;
                rep->neighbors(u, [&](int v, double peso) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        q.push(v);
                    }
                });
            }
            return pair<int,int>(ultimo, dist[ultimo]);
        };
        int inicio = 1;
        auto r1 = bfsMaisDistante(inicio);
        auto r2 = bfsMaisDistante(r1.first);
        return r2.second;
    }

    // Dijkstra com vetor simples (O(n²))
    void dijkstraVetor(int origem, const string& arquivo) const {
        const double INF = 1e18;
        vector<double> dist(n + 1, INF);
        vector<int> pai(n + 1, -1);
        vector<bool> visitado(n + 1, false);

        // Verifica se há peso negativo
        bool negativo = false;
        for (int u = 1; u <= n && !negativo; ++u) {
            rep->neighbors(u, [&](int v, double peso) {
                if (peso < 0) negativo = true;
            });
        }
        if (negativo) {
            cerr << "Erro: o grafo possui pesos negativos. Dijkstra não é aplicável.\n";
            return;
        }

        dist[origem] = 0.0;

        for (int i = 1; i <= n; ++i) {
            // Encontra vértice não visitado com menor distância
            int u = -1;
            double menor = INF;
            for (int v = 1; v <= n; ++v) {
                if (!visitado[v] && dist[v] < menor) {
                    menor = dist[v];
                    u = v;
                }
            }

            if (u == -1) break; // todos alcançáveis já processados
            visitado[u] = true;

            // Relaxa arestas de u
            rep->neighbors(u, [&](int v, double peso) {
                if (!visitado[v] && dist[u] + peso < dist[v]) {
                    dist[v] = dist[u] + peso;
                    pai[v] = u;
                }
            });
        }

        // Gera arquivo de saída
        //ofstream out(arquivo);
        //out << fixed << setprecision(2);
        //out << "Dijkstra (vetor) - Origem: " << origem << "\n";
        //out << "Vertice\tDistancia\tPai\tCaminho\n";
        //for (int i = 1; i <= n; ++i) {
            //if (dist[i] == INF) {
                //out << i << "\tInfinito\t" << pai[i] << "\tN/A\n";
                //continue;
            //}
            //out << i << "\t" << dist[i] << "\t\t" << pai[i] << "\t";

            // Reconstrói caminho
            //vector<int> caminho;
            //for (int v = i; v != -1; v = pai[v])
                //caminho.push_back(v);
            //reverse(caminho.begin(), caminho.end());
            //for (size_t j = 0; j < caminho.size(); ++j) {
                //out << caminho[j];
                //if (j + 1 < caminho.size()) out << " -> ";
            //}
            //out << "\n";
        //}
        //out.close();
    }

    void dijkstraHeap(int origem, const string& arquivo) const {
        vector<double> dist(n + 1, 1e18);
        vector<int> pai(n + 1, -1);
        vector<bool> visitado(n + 1, false);
        MinHeap heap(n);

        // inicializa
        dist[origem] = 0.0;
        heap.push(origem, 0.0);

        while (!heap.empty()) {
            auto minNode = heap.extract_min();
            int u = minNode.first;
            if (visitado[u]) continue;
            visitado[u] = true;

            rep->neighbors(u, [&](int v, double peso) {
                if (peso < 0) {
                    cerr << "Erro: peso negativo detectado entre " << u << " e " << v << endl;
                    throw runtime_error("Dijkstra não suporta pesos negativos.");
                }
                if (!visitado[v] && dist[u] + peso < dist[v]) {
                    dist[v] = dist[u] + peso;
                    pai[v] = u;
                    heap.push(v, dist[v]);
                }
            });
        }

        // salva resultados
        //ofstream out(arquivo);
        //out << fixed << setprecision(2);
        //out << "Dijkstra (heap) - Origem: " << origem << "\n";
        //out << "Vertice\tDistancia\tPai\tCaminho\n";

        //for (int i = 1; i <= n; ++i) {
            //out << i << "\t" << dist[i] << "\t\t" << pai[i] << "\t";
            //if (dist[i] == 1e18) { out << "(inacessível)\n"; continue; }

            // reconstrói caminho
            //vector<int> caminho;
            //for (int v = i; v != -1; v = pai[v]) caminho.push_back(v);
            //reverse(caminho.begin(), caminho.end());
            //for (size_t j = 0; j < caminho.size(); ++j) {
                //out << caminho[j];
                //if (j + 1 < caminho.size()) out << " -> ";
            //}
            //out << "\n";
        //}

        //out.close();
    }

    pair<vector<double>, vector<int>> dijkstraComResultado(int origem) const {
        const double INF = 1e18;
        vector<double> dist(n + 1, INF);
        vector<int> pai(n + 1, -1);
        vector<bool> visitado(n + 1, false);
        MinHeap heap(n);

        dist[origem] = 0.0;
        heap.push(origem, 0.0);

        while (!heap.empty()) {
            auto minNode = heap.extract_min();
            int u = minNode.first;
            if (visitado[u]) continue;
            visitado[u] = true;

            rep->neighbors(u, [&](int v, double peso) {
                if (!visitado[v] && dist[u] + peso < dist[v]) {
                    dist[v] = dist[u] + peso;
                    pai[v] = u;
                    heap.push(v, dist[v]);
                }
            });
        }

        return {dist, pai};
    }

    // Reconstrói o caminho do destino até a origem
    vector<int> reconstruirCaminho(int destino, const vector<int>& pai) const {
        vector<int> caminho;
        for (int v = destino; v != -1; v = pai[v]) {
            caminho.push_back(v);
        }
        reverse(caminho.begin(), caminho.end());
        return caminho;
    }
};

// Classe para gerenciar o mapeamento nome → vértice
class MapeamentoNomes {
private:
    map<string, int> nomeParaVertice;
    map<int, string> verticeParaNome;

public:
    // Lê o arquivo de mapeamento (formato: "id,Nome Completo")
    void lerMapeamento(const string& arquivo) {
        ifstream entrada(arquivo);
        if (!entrada.is_open()) {
            cerr << "Erro ao abrir arquivo de mapeamento: " << arquivo << endl;
            exit(1);
        }

        string linha;
        while (getline(entrada, linha)) {
            // Encontra a posição da vírgula
            size_t pos = linha.find(',');
            if (pos == string::npos) continue;

            // Extrai ID e nome
            int id = stoi(linha.substr(0, pos));
            string nome = linha.substr(pos + 1);

            // Remove espaços em branco do início e fim do nome
            nome.erase(0, nome.find_first_not_of(" \t\r\n"));
            nome.erase(nome.find_last_not_of(" \t\r\n") + 1);

            nomeParaVertice[nome] = id;
            verticeParaNome[id] = nome;
        }

        entrada.close();
        cout << "Mapeamento carregado: " << nomeParaVertice.size() << " pesquisadores." << endl;
    }

    int getVertice(const string& nome) const {
        auto it = nomeParaVertice.find(nome);
        if (it != nomeParaVertice.end()) {
            return it->second;
        }
        return -1; // Não encontrado
    }

    string getNome(int vertice) const {
        auto it = verticeParaNome.find(vertice);
        if (it != verticeParaNome.end()) {
            return it->second;
        }
        return "Desconhecido";
    }

    bool existe(const string& nome) const {
        return nomeParaVertice.find(nome) != nomeParaVertice.end();
    }
};

// ---------- main de exemplo ----------
int main() {
    // Lê o grafo do arquivo
    //Grafo g = Grafo::lerDeArquivo("grafo_W_5.txt", LISTA);

    //int origem = 10;

    //cout << "Executando Dijkstra (heap) a partir do vertice " << origem << "...\n";

    // Executa Dijkstra com heap e salva resultado em arquivo
    //g.dijkstraHeap(origem, "resultado_dijkstra_heap.txt");


    // Lê o grafo
    // Grafo g = Grafo::lerDeArquivo("grafo_W_2.txt", LISTA);
    // cout << "Grafo carregado."<< "\n";

    // // Quantos vértices aleatórios testar
    // int k = 100;
    // int totalVertices = 25000; // ajuste conforme o grafo real
    // vector<int> vertices;

    // // Gera k vértices aleatórios distintos entre 1 e totalVertices
    // mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    // uniform_int_distribution<int> dist(1, totalVertices);

    // for (int i = 0; i < k; ++i) {
    //     vertices.push_back(dist(rng));
    // }

    // double tempoVetorTotal = 0.0, tempoHeapTotal = 0.0;

    // for (int origem : vertices) {
    //     // Mede tempo para Dijkstra com vetor
    //     auto ini = std::chrono::high_resolution_clock::now();
    //     g.dijkstraVetor(origem, "");
    //     auto fim = std::chrono::high_resolution_clock::now();
    //     tempoVetorTotal += std::chrono::duration<double, std::milli>(fim - ini).count();

    //     // Mede tempo para Dijkstra com heap
    //     ini = std::chrono::high_resolution_clock::now();
    //     g.dijkstraHeap(origem, "");
    //     fim = std::chrono::high_resolution_clock::now();
    //     tempoHeapTotal += std::chrono::duration<double, std::milli>(fim - ini).count();
    // }

    // cout << fixed << setprecision(3);
    // cout << "\n==== RESULTADOS MÉDIOS ====\n";
    // cout << "Número de execuções (k): " << k << "\n";
    // cout << "Tempo médio (Dijkstra vetor): " << (tempoVetorTotal / k) << " ms\n";
    // cout << "Tempo médio (Dijkstra heap):  " << (tempoHeapTotal / k) << " ms\n";
    // cout << "Aceleração (vetor / heap):   " << (tempoVetorTotal / tempoHeapTotal) << "x\n";

    // Carrega o mapeamento de nomes
    MapeamentoNomes mapa;
    mapa.lerMapeamento("rede_colaboracao_vertices.txt");

    // Carrega o grafo
    Grafo g = Grafo::lerDeArquivo("rede_colaboracao.txt", LISTA);
    cout << "Grafo carregado.\n\n";

    // Define os pesquisadores
    string origem_nome = "Edsger W. Dijkstra";
    vector<string> destinos_nomes = {
        "Alan M. Turing",
        "J. B. Kruskal",
        "Jon M. Kleinberg",
        "Éva Tardos",
        "Daniel R. Figueiredo"
    };

    // Verifica se a origem existe
    if (!mapa.existe(origem_nome)) {
        cerr << "Erro: Pesquisador de origem '" << origem_nome << "' não encontrado!" << endl;
        return 1;
    }

    int origem = mapa.getVertice(origem_nome);
    cout << "Origem: " << origem_nome << " (vértice " << origem << ")\n\n";

    // Executa Dijkstra uma vez a partir da origem
    auto resultado = g.dijkstraComResultado(origem);
    const auto& distancias = resultado.first;
    const auto& pais = resultado.second;

    // Para cada destino, exibe distância e caminho
    cout << "===== RESULTADOS =====\n\n";
    
    for (const string& destino_nome : destinos_nomes) {
        cout << "Destino: " << destino_nome << "\n";
        
        if (!mapa.existe(destino_nome)) {
            cout << "  >> Pesquisador não encontrado no mapeamento.\n\n";
            continue;
        }

        int destino = mapa.getVertice(destino_nome);
        double dist = distancias[destino];

        if (dist >= 1e17) {
            cout << "  >> Não há caminho entre " << origem_nome << " e " << destino_nome << "\n\n";
            continue;
        }

        cout << "  Vértice: " << destino << "\n";
        cout << "  Distância: " << fixed << setprecision(2) << dist << "\n";
        cout << "  Caminho: ";

        vector<int> caminho = g.reconstruirCaminho(destino, pais);
        for (size_t i = 0; i < caminho.size(); ++i) {
            cout << mapa.getNome(caminho[i]);
            if (i + 1 < caminho.size()) cout << " -> ";
        }
        cout << "\n\n";
    }

    return 0;
}

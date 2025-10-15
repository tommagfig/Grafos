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
using namespace std;

enum Representacao { LISTA, MATRIZ };

struct IGrafoRep {
    virtual ~IGrafoRep() = default;
    virtual void add_edge(int u, int v) = 0;
    virtual void neighbors(int u, const function<void(int)>& f) const = 0;
    virtual int degree(int u) const = 0;
    virtual bool edge_exists(int u, int v) const = 0;
    virtual int n_vertices() const = 0;
};

// Nó da lista de adjacência
struct Node {
    int v;
    Node* next;
    Node(int vert, Node* nxt = nullptr) : v(vert), next(nxt) {}
};

//lista de adjacência
class ListaRep : public IGrafoRep {
private:
    vector<Node*> adj; 

public:
    explicit ListaRep(int n) : adj(n + 1, nullptr) {}

    void add_edge(int u, int v) override {
        adj[u] = new Node(v, adj[u]);
        adj[v] = new Node(u, adj[v]);
    }

    void neighbors(int u, const function<void(int)>& f) const override {
        for (Node* curr = adj[u]; curr != nullptr; curr = curr->next) {
            f(curr->v);
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
    vector<vector<uint8_t>> mat;
public:
    explicit MatrizRep(int n) : mat(n + 1, vector<uint8_t>(n + 1, 0)) {}
    void add_edge(int u, int v) override {
        mat[u][v] = 1;
        mat[v][u] = 1;
    }
    void neighbors(int u, const function<void(int)>& f) const override {
        int n = (int)mat[u].size() - 1;
        for (int v = 1; v <= n; ++v) {
            if (mat[u][v]) f(v);
        }
    }
    int degree(int u) const override {
        int cnt = 0;
        int n = (int)mat[u].size() - 1;
        for (int v = 1; v <= n; ++v) if (mat[u][v]) ++cnt;
        return cnt;
    }
    bool edge_exists(int u, int v) const override {
        return mat[u][v] != 0;
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

    void adicionarAresta(int u, int v) {
        if (u < 1 || u > n || v < 1 || v > n) return;
        if (!rep->edge_exists(u,v)) {
            rep->add_edge(u, v);
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
        entrada >> vertices;
        Grafo g(vertices, r);
        while (entrada >> u >> v) {
            g.adicionarAresta(u, v);
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
            rep->neighbors(u, [&](int v){
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
                rep->neighbors(u, [&](int v){
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
                    rep->neighbors(u, [&](int v){
                        if (!visitado[v]) { visitado[v] = 1; q.push(v); }
                    });
                }
                componentes.push_back(move(comp));
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
            rep->neighbors(u, [&](int v){
                if (nivel[v] == -1) {
                    pai[v] = u;
                    nivel[v] = nivel[u] + 1;
                    q.push(v);
                }
            });
        }

        ofstream out(arquivo);
        out << "BFS - Arvore de busca a partir do vertice " << origem << "\n";
        out << "Vertice\tPai\tNivel\n";
        for (int i = 1; i <= n; ++i) out << i << "\t\t" << pai[i] << "\t" << nivel[i] << "\n";
        out.close();
    }

    // DFS iterativa que escreve pai e nivel em arquivo
    void DFS(int origem, const string& arquivo) const {
        vector<int> nivel(n + 1, -1);
        vector<int> pai(n + 1, 0);
        stack<pair<int,int>> s;
        s.push({origem, 0});
        pai[origem] = 0;

        while (!s.empty()) {
            auto [u, d] = s.top(); s.pop();
            if (nivel[u] != -1) continue;
            nivel[u] = d;

            vector<int> vizinhos;
            rep->neighbors(u, [&](int v){ vizinhos.push_back(v); });
            for (auto it = vizinhos.rbegin(); it != vizinhos.rend(); ++it) {
                int v = *it;
                if (nivel[v] == -1) {
                    pai[v] = u;
                    s.push({v, d + 1});
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
                rep->neighbors(u, [&](int v){
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

    
};

// ---------- main de exemplo (idêntico ao anterior) ----------
int main() {
    
}

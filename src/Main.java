public class Main {


    public static void main(String[] args){

        Vertex v1 = new Vertex("Vancouver",80,78);
        Vertex v2 = new Vertex("Calgary",180,67);
        Vertex v3 = new Vertex("Winnipeg",360,75);
        Vertex v4 = new Vertex("Montreal",696,66);
        Vertex v5 = new Vertex("SaultSaintMarie",550,120);

        Dijkstra bleh = new Dijkstra();

        bleh.addVertex(v1);
        bleh.addVertex(v2);
        bleh.addVertex(v3);
        bleh.addVertex(v4);
        bleh.addVertex(v5);

        bleh.computeAllEuclideanDistances();

    }

}
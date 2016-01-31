#include "DigraphExample.h"

/**
 * Used to test deep first search
 * The digraph comes from Figure 22-4 of "Cormen et al: Introduction to
 * agorithms", Chapter 22.3.
 */
Digraph* DigraphExample::generateDigraph0() {
	Digraph* digraph = new Digraph();

	// creat Node
	Node* u = new Node("u", 1, 1,1,4);
	Node* v = new Node("v", 2, 1,2,4);
	Node* w = new Node("w", 3, 1,2,4);
	Node* x = new Node("x", 4, 1,3,4);
	Node* y = new Node("y", 5, 1,2,4);
	Node* z = new Node("z", 6, 1,3,4);

	digraph->add_node(u);
	digraph->add_node(v);
	digraph->add_node(w);
	digraph->add_node(x);
	digraph->add_node(y);
	digraph->add_node(z);

	// Create Edge
	Edge* uv = new Edge(u,v);
	Edge* ux = new Edge(u,x);
	Edge* vy = new Edge(v,y);
	Edge* wy = new Edge(w,y);
	Edge* wz = new Edge(w,z);
	Edge* xv = new Edge(x,v);
	Edge* yx = new Edge(y,x);
	Edge* zz = new Edge(z,z);

	uv->set_separation_time(4);
	ux->set_separation_time(4);
	vy->set_separation_time(4);
	wy->set_separation_time(4);
	wz->set_separation_time(4);
	xv->set_separation_time(4);
	yx->set_separation_time(4);
	zz->set_separation_time(4);

	digraph->add_edge(uv);
	digraph->add_edge(ux);
	digraph->add_edge(vy);
	digraph->add_edge(wy);
	digraph->add_edge(wz);
	digraph->add_edge(xv);
	digraph->add_edge(yx);
	digraph->add_edge(zz);
	
	if (false) {
		digraph->write_graphviz(std::cout);
	}

	return digraph;
}


/**
 * Used to test strongly connected components algortihm
 * The digraph comes from Figure 22-9 of "Cormen et al: Introduction to
 * agorithms", Chapter 22.5. 
 */
Digraph* DigraphExample::generateDigraph1() {
	Digraph* digraph = new Digraph();

	// creat Node
	Node* a = new Node("a", 1, 10);
	Node* b = new Node("b", 2, 10);
	Node* c = new Node("c", 3, 10);
	Node* d = new Node("d", 4, 10);
	Node* e = new Node("e", 5, 10);
	Node* f = new Node("f", 6, 10);
	Node* g = new Node("g", 7, 10);
	Node* h = new Node("h", 8, 10);

	digraph->add_node(a);
	digraph->add_node(b);
	digraph->add_node(c);
	digraph->add_node(d);
	digraph->add_node(e);
	digraph->add_node(f);
	digraph->add_node(g);
	digraph->add_node(h);

	// Create Edge
	Edge* ab = new Edge(a,b);
	Edge* bf = new Edge(b,f);
	Edge* be = new Edge(b,e);
	Edge* ea = new Edge(e,a);
	Edge* ef = new Edge(e,f);
	Edge* bc = new Edge(b,c);
	Edge* fg = new Edge(f,g);
	Edge* gf = new Edge(g,f);
	Edge* cd = new Edge(c,d);
	Edge* dc = new Edge(d,c);
	Edge* cg = new Edge(c,g);
	Edge* gh = new Edge(g,h);
	Edge* hh = new Edge(h,h);
	Edge* dh = new Edge(d,h);

	digraph->add_edge(ab);
	digraph->add_edge(bf);
	digraph->add_edge(be);
	digraph->add_edge(ea);
	digraph->add_edge(ef);
	digraph->add_edge(bc);
	digraph->add_edge(fg);
	digraph->add_edge(gf);
	digraph->add_edge(cd);
	digraph->add_edge(dc);
	digraph->add_edge(cg);
	digraph->add_edge(gh);
	digraph->add_edge(hh);
	digraph->add_edge(dh);

	if (false) {
		digraph->write_graphviz(std::cout);
	}

	return digraph;
}

/**
 * An example of the paper RTNS2015.
 * Used to calculate the linear factor (or utilization)
 */
Digraph* DigraphExample::generateDigraph2() {
	Digraph* digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,2,3);
	Node* v2 = new Node("v2",1,2,3);
	Node* v3 = new Node("v3",1,1,2);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);

	// create edges
	Edge* v1_v1 = new Edge(v1,v1);
	v1_v1->set_separation_time(3);

	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(4);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(3);

	Edge* v3_v2 = new Edge(v3,v2);
	v3_v2->set_separation_time(2);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(2);

	digraph->add_edge(v1_v1);
	digraph->add_edge(v1_v2);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v3_v2);
	digraph->add_edge(v3_v1);

	if (false) {
		digraph->write_graphviz(std::cout);
	}

	return digraph;
}

/**
 * An example of the paper RTS2014.
 * Used to test rbf and dbf
 */
Digraph* DigraphExample::generateDigraph3() {
	Digraph* digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,1,10);
	Node* v2 = new Node("v2",1,2,10);
	Node* v3 = new Node("v3",1,1,10);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);

	// create edges
	Edge* v1_v1 = new Edge(v1,v1);
	v1_v1->set_separation_time(10);

	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(20);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(10);

	Edge* v3_v2 = new Edge(v3,v2);
	v3_v2->set_separation_time(20);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(10);

	digraph->add_edge(v1_v1);
	digraph->add_edge(v1_v2);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v3_v2);
	digraph->add_edge(v3_v1);

	if (false) {
		digraph->write_graphviz(std::cout);
	}

	return digraph;
}

/**
 * An example of Martin Stigge et al., The digraph real-time task model, RTAS2011
 * Used to dbf
 */
Digraph* DigraphExample::generateDigraph4() {
	Digraph* digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,2,5);
	Node* v2 = new Node("v2",1,1,8);
	Node* v3 = new Node("v3",1,3,8);
	Node* v4 = new Node("v4",1,5,10);
	Node* v5 = new Node("v5",1,1,5);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);
	digraph->add_node(v4);
	digraph->add_node(v5);

	// create edges
	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(10);

	Edge* v1_v5 = new Edge(v1,v5);
	v1_v5->set_separation_time(20);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(15);

	Edge* v2_v4 = new Edge(v2,v4);
	v2_v4->set_separation_time(20);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(11);

	Edge* v4_v2 = new Edge(v4,v2);
	v4_v2->set_separation_time(20);

	Edge* v5_v4 = new Edge(v5,v4);
	v5_v4->set_separation_time(10);

	digraph->add_edge(v1_v2);
	digraph->add_edge(v1_v5);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v2_v4);
	digraph->add_edge(v3_v1);
	digraph->add_edge(v4_v2);
	digraph->add_edge(v5_v4);

	if (false) {
		digraph->write_graphviz(std::cout);
	}

	return digraph;
}

Digraph* DigraphExample::generateDigraph5() {
	Digraph* digraph = new Digraph();

	// create nodes
	Node* A10 = new Node("A10",1,10,200);
	Node* A12 = new Node("A12",1,10,200);
	Node* A14 = new Node("A14",1,10,100);
	Node* A16 = new Node("A16",1,10,200);
	Node* A18 = new Node("A18",1,10,200);

	Node* A30 = new Node("A30",1,25,500);
	Node* A32 = new Node("A32",1,25,300);
	Node* A34 = new Node("A34",1,25,100);
	Node* A36 = new Node("A36",1,25,400);
	Node* A38 = new Node("A38",1,25,200);

	Node* A40 = new Node("A40",1,15,500);
	Node* A45 = new Node("A45",1,15,500);

	Node* A20 = new Node("A20",1,30,200);
	Node* A25 = new Node("A25",1,30,100);

	digraph->add_node(A10);
	digraph->add_node(A12);
	digraph->add_node(A14);
	digraph->add_node(A16);
	digraph->add_node(A18);

	digraph->add_node(A30);
	digraph->add_node(A32);
	digraph->add_node(A34);
	digraph->add_node(A36);
	digraph->add_node(A38);

	digraph->add_node(A40);
	digraph->add_node(A45);

	digraph->add_node(A20);
	digraph->add_node(A25);

	// create edges
	Edge* A10_A32 = new Edge(A10,A32);
	A10_A32->set_separation_time(200);

	Edge* A10_A45 = new Edge(A10,A45);
	A10_A45->set_separation_time(500);

	Edge* A12_A34 = new Edge(A12,A34);
	A12_A34->set_separation_time(200);

	Edge* A12_A45 = new Edge(A12,A45);
	A12_A45->set_separation_time(300);

	Edge* A14_A45 = new Edge(A14,A45);
	A14_A45->set_separation_time(100);

	Edge* A14_A36 = new Edge(A14,A36);
	A14_A36->set_separation_time(200);

	Edge* A16_A38 = new Edge(A16,A38);
	A16_A38->set_separation_time(200);

	Edge* A16_A40 = new Edge(A16,A40);
	A16_A40->set_separation_time(400);

	Edge* A18_A30 = new Edge(A18,A30);
	A18_A30->set_separation_time(200);

	Edge* A18_A40 = new Edge(A18,A40);
	A18_A40->set_separation_time(200);

	Edge* A30_A25 = new Edge(A30,A25);
	A30_A25->set_separation_time(500);

	Edge* A32_A25 = new Edge(A32,A25);
	A32_A25->set_separation_time(300);

	Edge* A34_A25 = new Edge(A34,A25);
	A34_A25->set_separation_time(100);

	Edge* A36_A20 = new Edge(A36,A20);
	A36_A20->set_separation_time(400);

	Edge* A38_A20 = new Edge(A38,A20);
	A38_A20->set_separation_time(200);

	Edge* A40_A25 = new Edge(A40,A25);
	A40_A25->set_separation_time(500);

	Edge* A45_A20 = new Edge(A45,A20);
	A45_A20->set_separation_time(500);

	Edge* A20_A12 = new Edge(A20,A12);
	A20_A12->set_separation_time(200);

	Edge* A25_A16 = new Edge(A25,A16);
	A25_A16->set_separation_time(100);

	digraph->add_edge(A10_A32);
	digraph->add_edge(A10_A45);
	digraph->add_edge(A12_A34);
	digraph->add_edge(A12_A45);
	digraph->add_edge(A14_A45);
	digraph->add_edge(A14_A36);
	digraph->add_edge(A16_A38);
	digraph->add_edge(A16_A40);
	digraph->add_edge(A18_A30);
	digraph->add_edge(A18_A40);
	digraph->add_edge(A30_A25);
	digraph->add_edge(A32_A25);
	digraph->add_edge(A34_A25);
	digraph->add_edge(A36_A20);
	digraph->add_edge(A38_A20);
	digraph->add_edge(A40_A25);
	digraph->add_edge(A45_A20);
	digraph->add_edge(A20_A12);
	digraph->add_edge(A25_A16);

	return digraph;
}

Digraph* DigraphExample::generateDigraph6() {
	Digraph* digraph = new Digraph();

	// create nodes
	Node* A10 = new Node("A10",1,10,200);
	Node* A12 = new Node("A12",1,10,200);
	Node* A14 = new Node("A14",1,10,100);
	Node* A16 = new Node("A16",1,10,200);
	Node* A18 = new Node("A18",1,10,200);

	Node* A30 = new Node("A30",1,25,500);
	Node* A32 = new Node("A32",1,25,300);
	Node* A34 = new Node("A34",1,25,100);
	Node* A36 = new Node("A36",1,25,400);
	Node* A38 = new Node("A38",1,25,200);

	Node* A40 = new Node("A40",1,15,500);
	Node* A45 = new Node("A45",1,15,500);

	Node* A20 = new Node("A20",1,30,200);
	Node* A25 = new Node("A25",1,30,100);

	digraph->add_node(A10);
	digraph->add_node(A12);
	digraph->add_node(A14);
	digraph->add_node(A16);
	digraph->add_node(A18);

	digraph->add_node(A30);
	digraph->add_node(A32);
	digraph->add_node(A34);
	digraph->add_node(A36);
	digraph->add_node(A38);

	digraph->add_node(A40);
	digraph->add_node(A45);

	digraph->add_node(A20);
	digraph->add_node(A25);

	// create edges
	Edge* A10_A32 = new Edge(A10,A32);
	A10_A32->set_separation_time(200);

	Edge* A10_A45 = new Edge(A10,A45);
	A10_A45->set_separation_time(500);

	Edge* A12_A34 = new Edge(A12,A34);
	A12_A34->set_separation_time(200);

	Edge* A12_A45 = new Edge(A12,A45);
	A12_A45->set_separation_time(300);

	Edge* A14_A45 = new Edge(A14,A45);
	A14_A45->set_separation_time(100);

	Edge* A14_A36 = new Edge(A14,A36);
	A14_A36->set_separation_time(200);

	Edge* A16_A38 = new Edge(A16,A38);
	A16_A38->set_separation_time(200);

	Edge* A16_A40 = new Edge(A16,A40);
	A16_A40->set_separation_time(400);

	Edge* A18_A30 = new Edge(A18,A30);
	A18_A30->set_separation_time(200);

	Edge* A18_A40 = new Edge(A18,A40);
	A18_A40->set_separation_time(200);

	Edge* A30_A25 = new Edge(A30,A25);
	A30_A25->set_separation_time(500);

	Edge* A32_A25 = new Edge(A32,A25);
	A32_A25->set_separation_time(300);

	Edge* A34_A25 = new Edge(A34,A25);
	A34_A25->set_separation_time(100);

	Edge* A36_A20 = new Edge(A36,A20);
	A36_A20->set_separation_time(400);

	Edge* A38_A20 = new Edge(A38,A20);
	A38_A20->set_separation_time(200);

	Edge* A40_A25 = new Edge(A40,A25);
	A40_A25->set_separation_time(500);

	Edge* A45_A20 = new Edge(A45,A20);
	A45_A20->set_separation_time(500);

	Edge* A20_A12 = new Edge(A20,A12);
	A20_A12->set_separation_time(200);

	Edge* A25_A16 = new Edge(A25,A16);
	A25_A16->set_separation_time(100);

	digraph->add_edge(A10_A32);
	digraph->add_edge(A10_A45);
	digraph->add_edge(A12_A34);
	digraph->add_edge(A12_A45);
	digraph->add_edge(A14_A45);
	digraph->add_edge(A14_A36);
	digraph->add_edge(A16_A38);
	digraph->add_edge(A16_A40);
	digraph->add_edge(A18_A30);
	digraph->add_edge(A18_A40);
	digraph->add_edge(A30_A25);
	digraph->add_edge(A32_A25);
	digraph->add_edge(A34_A25);
	digraph->add_edge(A36_A20);
	digraph->add_edge(A38_A20);
	digraph->add_edge(A40_A25);
	digraph->add_edge(A45_A20);
	digraph->add_edge(A20_A12);
	digraph->add_edge(A25_A16);

	// add some loose edges such that the modified digraph is strongly connected
	Edge* A25_A18 = new Edge(A25,A18);
	A25_A18->set_separation_time(300);

	Edge* A25_A10 = new Edge(A25,A10);
	A25_A10->set_separation_time(500);

	Edge* A20_A14 = new Edge(A20,A14);
	A20_A14->set_separation_time(400);

	digraph->add_edge(A25_A18);
	digraph->add_edge(A25_A10);
	digraph->add_edge(A20_A14);

	return digraph;
}

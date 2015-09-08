#include "NodeAndEdge.h"

// ===========================================================================
// Method functions of Node class
// ===========================================================================
Node::Node(std::string _name) {
	name = _name;
}

Node::Node(std::string _name, int _index, int _scale) {
	name = _name;
	index = _index;
	scale = _scale;
}

Node::Node(std::string _name, int _scale, int _wcet, int _deadline) {
	name = _name;
	scale = _scale;

	set_wcet(_wcet);
	set_deadline(_deadline);
}

Node::Node(std::string _name, int _index, int _scale, int _wcet, int _deadline) {
	name = _name;
	index = _index;
	scale = _scale;

	set_wcet(_wcet);
	set_deadline(_deadline);
}

void Node::set_wcet(int _wcet) {
	wcet = _wcet;
}

void Node::set_deadline(int _deadline) {
	deadline = _deadline;
}

std::string Node::toString() {
	return name;
}

void Node::set_color(Color _color) {
	this->color = _color;
}

Node::Color Node::get_color() {
	return color;
}
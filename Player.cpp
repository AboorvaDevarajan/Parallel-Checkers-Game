#include "Player.h"

Player::Player() {
	name = "";
	isMax = true;
}
Player::Player(string myname, bool maxMin) {
	name = myname;
	isMax = maxMin;
}

Player::~Player() {

}


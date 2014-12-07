#include <vector>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <vector>
#include <time.h>
#include <limits>



#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif


#include <Eigen/Dense>
#include <Eigen/StdVector>

using namespace std;
using namespace Eigen;


const double THETA_CHANGE = 0.1;


Vector3d rotateX(const Vector3d& p, double theta) {
    Matrix3d rotateX;
    rotateX <<
        1, 0, 0,
        0, cos(theta), -sin(theta),
        0, sin(theta), cos(theta);

    return rotateX * p;
}

Vector3d rotateY(const Vector3d& p, double theta) {
    Matrix3d rotateY;
    rotateY <<
        cos(theta), 0, sin(theta),
        0, 1, 0,
        -sin(theta), 0, cos(theta);

    return rotateY * p;
}

Vector3d rotateZ(const Vector3d& p, double theta) {
    Matrix3d rotateZ;
    rotateZ <<
        cos(theta), -sin(theta), 0,
        sin(theta), cos(theta), 0,
        0, 0, 1;

    return rotateZ * p;
}


Vector3d translate(const Vector3d& p, const Vector3d& to) {
    return p - to;
}



class Viewport {
public:
    int w, h; // width and height
};


class Joint {
public:
    double length;
    Vector3d theta;
    Vector3d basePoint, endPoint;

    void update(Vector3d dtheta, Vector3d newBasePoint);
    void render();

    Joint(Vector3d inBasePoint, double inLength);
};


Joint::Joint(Vector3d inBasePoint, double inLength) {
    length = inLength;
    theta = Vector3d(0.0, 0.0, 0.0);
    basePoint = inBasePoint;
    endPoint = basePoint + Vector3d(0.0, length, 0.0);
}

void Joint::update(Vector3d dtheta, Vector3d newBasePoint) {
    basePoint = newBasePoint;

    Vector3d original(0.0, length, 0.0);

    theta += dtheta;

    Vector3d newEndPoint = rotateZ(rotateY(rotateX(original, theta[0]), theta[1]), theta[2]);
    endPoint = translate(newEndPoint, -newBasePoint);

}

void Joint::render() {
    glBegin(GL_LINES);
        glVertex3d(basePoint[0], basePoint[1], basePoint[2]);
        glVertex3d(endPoint[0], endPoint[1], endPoint[2]);
    glEnd();
}



class System {
private:
    vector<Joint> joints;
    MatrixXd getJacobian();
    void updateJoints(const VectorXd& dtheta);
public:

    System(vector<Joint> inJoints);
    bool update(Vector3d g);
    void render();

};


System::System(vector<Joint> inJoints) {
    joints = inJoints;
}


MatrixXd System::getJacobian() {
    MatrixXd result(3, 3 * joints.size());

    Vector3d deltaTheta(1.0, 1.0, 1.0);
    deltaTheta *= THETA_CHANGE;

    for (int i = 0; i < joints.size(); i++) {
        Vector3d newTheta = joints[i].theta + deltaTheta;

        Vector3d original(0.0, joints[i].length, 0.0); //translate(joints[i].endPoint, joints[i].basePoint);


        Vector3d p = rotateX(original, joints[i].theta[0]);
        Vector3d newP = rotateX(original, newTheta[0]);

        Vector3d col1 = (newP - p) / THETA_CHANGE;

        result(0, i * 3) = col1[0];
        result(1, i * 3) = col1[1];
        result(2, i * 3) = col1[2];

        p = rotateY(original, joints[i].theta[1]);
        newP = rotateY(original, newTheta[1]);

        Vector3d col2 = (newP - p) / THETA_CHANGE;

        result(0, i * 3 + 1) = col2[0];
        result(1, i * 3 + 1) = col2[1];
        result(2, i * 3 + 1) = col2[2];

        p = rotateZ(original, joints[i].theta[2]);
        newP = rotateZ(original, newTheta[2]);

        Vector3d col3 = (newP - p) / THETA_CHANGE;

        result(0, i * 3 + 2) = col3[0];
        result(1, i * 3 + 2) = col3[1];
        result(2, i * 3 + 2) = col3[2];
    }


    return result;
}


void System::updateJoints(const VectorXd& dtheta) {

    Vector3d newBasePoint = joints[0].basePoint;
    Vector3d currDTheta(0.0, 0.0, 0.0);
    for(int i = 0; i < joints.size(); i++) {
        currDTheta = Vector3d(dtheta[3*i], dtheta[3*i+1], dtheta[3*i+2]);
        joints[i].update(currDTheta, newBasePoint);
        newBasePoint = joints[i].endPoint;
    }
}


bool System::update(Vector3d g) {
    Vector3d gSys = g - joints[0].basePoint;

    double jointsLength = 0.0;

    for (int i = 0; i < joints.size(); i++) {
        jointsLength += joints[i].length;
    }

    if (gSys.norm() > jointsLength) {
        cout << "not reachable" << endl;
    }
    Vector3d dp = g - joints[joints.size() - 1].endPoint;
    if (dp.norm() > 0.0001) {
        MatrixXd J = getJacobian();

        VectorXd dtheta = J.jacobiSvd(ComputeThinU | ComputeThinV).solve(dp);

        updateJoints(dtheta);
        return true;
    }

    return false;
}


double colors[] = {0.0, 0.0, 1.0,
                   0.0, 1.0, 0.0,
                   1.0, 0.0, 0.0 };


void System::render() {
    for(int i = 0; i < joints.size(); i++) {
        glColor3d(colors[i * 3], colors[i * 3 + 1], colors[i * 3 + 2]);
        joints[i].render();
        // cout << (joints[i].endPoint - joints[i].basePoint).norm() << "lol" << endl;
        cout << joints[i].endPoint << "lol" << endl;
    }
}



Viewport viewport;
vector<Joint> bones;
System arm(bones);


void myReshape(int w, int h) {
    viewport.w = w;
    viewport.h = h;

    glViewport(0, 0, viewport.w, viewport.h);// sets the rectangle that will be the window
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();                // loading the identity matrix for the screen

    //----------- setting the projection -------------------------
    // glOrtho sets left, right, bottom, top, zNear, zFar of the chord system
    // glOrtho(-1, 1 + (w-400)/200.0 , -1 -(h-400)/200.0, 1, 1, -1); // resize type = add
     glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 2, -2); // resize type = center

    //glOrtho(-4, 4, -4, 4, 2, -2);    // resize type = stretch

    //------------------------------------------------------------
}


int lol;


void initScene(){
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent


    bones.push_back(Joint(Vector3d(0.0, 0.0, 0.0), 0.5));
    bones.push_back(Joint(bones[0].endPoint, 0.5));
    bones.push_back(Joint(bones[1].endPoint, 0.5));

    arm = System(bones);

    myReshape(viewport.w, viewport.h);
}



void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    glShadeModel(GL_FLAT);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Vector3d goal(0.7, 0.0, 0.7);


    arm.update(goal);
    arm.render();

    glFlush();
    glutSwapBuffers();
}

//****************************************************
// called by glut when there are no messages to handle
//****************************************************
void myFrameMove() {
    //nothing here for now
#ifdef _WIN32
    Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
#endif
    glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}


int main(int argc, char *argv[]) {

    glutInit(&argc, argv);

    //This tells glut to use a double-buffered window with red, green, and blue channels
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);


    // Initalize theviewport size
    viewport.w = 800;
    viewport.h = 800;

    //The size and position of the window
    glutInitWindowSize(viewport.w, viewport.h);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("CS184 Assignment 4");

    initScene();                                 // quick function to set up scene

    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);                    // function to run when its time to draw something
    glutReshapeFunc(myReshape);                  // function to run when the window gets resized
    glutIdleFunc(myFrameMove);                   // function to run when not handling any other task
    glutMainLoop();                              // infinite loop that will keep drawing and resizing and whatever else

    return 0;
}
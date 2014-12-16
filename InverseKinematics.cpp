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


const double SCALE_DELTA = 0.5;

const double THETA_CHANGE = 0.01;

const double PI = 3.141592654;


float zoomFactor = 1.0f;


Vector3d camera(0.0, 0.0, 10.0);

Matrix3d rotateX(double theta) {
    Matrix3d rotateX;
    rotateX <<
        1, 0, 0,
        0, cos(theta), -sin(theta),
        0, sin(theta), cos(theta);

    return rotateX;
}

Matrix3d rotateY(double theta) {
    Matrix3d rotateY;
    rotateY <<
        cos(theta), 0, sin(theta),
        0, 1, 0,
        -sin(theta), 0, cos(theta);

    return rotateY;
}


Matrix3d rotateZ(double theta) {
    Matrix3d rotateZ;
    rotateZ <<
        cos(theta), -sin(theta), 0,
        sin(theta), cos(theta), 0,
        0, 0, 1;

    return rotateZ;
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

    void update(Vector3d dtheta, Vector3d newBasePoint, Vector3d newEndPoint);
    void render(Vector3d thetaPrev);

    Joint(Vector3d inBasePoint, double inLength);
};


Joint::Joint(Vector3d inBasePoint, double inLength) {
    length = inLength;
    theta = Vector3d(0.0, 0.0, 0.0);
    basePoint = inBasePoint;
    endPoint = basePoint + Vector3d(0.0, length, 0.0);
}


void Joint::update(Vector3d dtheta, Vector3d newBasePoint, Vector3d newEndPoint) {
    basePoint = newBasePoint;
    theta += dtheta;
    endPoint = newEndPoint;
}


void Joint::render(Vector3d thetaPrev) {


    Vector3d base(0.0, 0.0, 0.0);
    Vector3d end(0.0, length, 0.0);


    Vector3d toDegrees = 180.0/PI * (thetaPrev + theta);

    double d = length / 6;
    Vector3d first    = base + Vector3d(-d, length/3, -d);
    Vector3d second   = base + Vector3d(d, length/3, -d);
    Vector3d third    = base + Vector3d(d, length/3, d);
    Vector3d fourth   = base + Vector3d(-d, length/3, d);

    glPushMatrix();

        glTranslated(basePoint[0], basePoint[1], basePoint[2]);
        glRotated(toDegrees[2], 0, 0, 1);
        glRotated(toDegrees[1], 0, 1, 0);
        glRotated(toDegrees[0], 1, 0, 0);

        glBegin(GL_POLYGON);
            glVertex3d(base[0], base[1], base[2]);
            glVertex3d(first[0], first[1], first[2]);
            glVertex3d(second[0], second[1], second[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(base[0], base[1], base[2]);
            glVertex3d(second[0], second[1], second[2]);
            glVertex3d(third[0], third[1], third[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(base[0], base[1], base[2]);
            glVertex3d(third[0], third[1], third[2]);
            glVertex3d(fourth[0], fourth[1], fourth[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(base[0], base[1], base[2]);
            glVertex3d(fourth[0], fourth[1], fourth[2]);
            glVertex3d(first[0], first[1], first[2]);
        glEnd();

        //second
        glBegin(GL_POLYGON);
            glVertex3d(end[0], end[1], end[2]);
            glVertex3d(second[0], second[1], second[2]);
            glVertex3d(first[0], first[1], first[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(end[0], end[1], end[2]);
            glVertex3d(third[0], third[1], third[2]);
            glVertex3d(second[0], second[1], second[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(end[0], end[1], end[2]);
            glVertex3d(fourth[0], fourth[1], fourth[2]);
            glVertex3d(third[0], third[1], third[2]);
        glEnd();

        glBegin(GL_POLYGON);
            glVertex3d(end[0], end[1], end[2]);
            glVertex3d(first[0], first[1], first[2]);
            glVertex3d(fourth[0], fourth[1], fourth[2]);
        glEnd();

    glPopMatrix();
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


    Vector3d p = joints[joints.size() - 1].endPoint;
    for (int i = 0; i < joints.size(); i++) {



        Vector3d original = translate(p, joints[i].basePoint);

        Vector3d newP = rotateX(THETA_CHANGE) * original;

        Vector3d col1 = (newP - original) / THETA_CHANGE;

        result(0, i * 3) = col1[0];
        result(1, i * 3) = col1[1];
        result(2, i * 3) = col1[2];

        newP = rotateY(THETA_CHANGE) * original;

        Vector3d col2 = (newP - original) / THETA_CHANGE;

        result(0, i * 3 + 1) = col2[0];
        result(1, i * 3 + 1) = col2[1];
        result(2, i * 3 + 1) = col2[2];

        newP = rotateZ(THETA_CHANGE) * original;

        Vector3d col3 = (newP - original) / THETA_CHANGE;

        result(0, i * 3 + 2) = col3[0];
        result(1, i * 3 + 2) = col3[1];
        result(2, i * 3 + 2) = col3[2];


    }

    return result;
}


void System::updateJoints(const VectorXd& dtheta) {

    Vector3d newBasePoint = joints[0].basePoint;
    Vector3d currDTheta(0.0, 0.0, 0.0);
    Vector3d currTheta(0.0, 0.0, 0.0);
    Vector3d original(0.0, 1.0, 0.0);


    for(int i = 0; i < joints.size(); i++) {
        currDTheta = Vector3d(dtheta[3*i], dtheta[3*i+1], dtheta[3*i+2]);
        currTheta += joints[i].theta + currDTheta;

        original = Vector3d(0.0, joints[i].length, 0.0);

        Vector3d newEndPoint = rotateZ(currTheta[2]) * rotateY(currTheta[1]) * rotateX(currTheta[0]) * original;
        newEndPoint = translate(newEndPoint, -newBasePoint);
        joints[i].update(currDTheta, newBasePoint, newEndPoint);
        newBasePoint = joints[i].endPoint;
    }
}



MatrixXd getPseudoInverse(MatrixXd& J) {
    MatrixXd transposeJ = J.transpose();

    return transposeJ * (J * transposeJ).inverse();
}



bool System::update(Vector3d g) {

    double jointsLength = 0.0;

    for (int i = 0; i < joints.size(); i++) {
        jointsLength += joints[i].length;
    }

    if (g.norm() > jointsLength) {
        g = g.normalized() * jointsLength;
    }

    Vector3d dp = g - joints[joints.size() - 1].endPoint;

    MatrixXd J = getJacobian();
    MatrixXd inverseJ = getPseudoInverse(J);

    VectorXd dtheta = J.jacobiSvd(ComputeThinU | ComputeThinV).solve(dp); //inverseJ * dp; //

    updateJoints(dtheta);

    return true;
}


double colors[] = {0.0, 0.0, 1.0,
                   0.0, 1.0, 0.0,
                   1.0, 0.0, 0.0,
                   1.0, 1.0, 0.0 };


void System::render() {
    Vector3d thetaPrev(0.0, 0.0, 0.0);
    for(int i = 0; i < joints.size(); i++) {
        GLfloat lightColor0[] = { colors[i * 3], colors[i * 3 + 1], colors[i * 3 + 2], 1.0f }; //Color (0.5, 0.5, 0.5)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
        joints[i].render(thetaPrev);
        thetaPrev += joints[i].theta;
    }
}


vector<Vector3d> goals;

int currGoalIndex;

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
    glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 2, -2); // resize type = center


    //------------------------------------------------------------
}


void processSpecialKeys(int key, int xx, int yy) {

    float fraction = 0.1f;


    int mod = glutGetModifiers();
    bool isShift = false;

    switch (mod) {
        case GLUT_ACTIVE_SHIFT:
            isShift = true;
            break;
    }

    switch (key) {
        case GLUT_KEY_LEFT :
            camera = rotateY(-fraction) * camera;
            break;
        case GLUT_KEY_RIGHT :
            camera = rotateY(fraction) * camera;
            break;
        case GLUT_KEY_UP :
            camera = rotateX(-fraction) * camera;
            break;
        case GLUT_KEY_DOWN :
            camera = rotateX(fraction) * camera;
            break;
    }

}



void myKeyboard (unsigned char key, int x, int y) {
    switch (key) {
        case '+':
            zoomFactor *= SCALE_DELTA;
            break;

        case '-':
            zoomFactor /= SCALE_DELTA;
            break;
    }

    glutPostRedisplay();
}



clock_t begin1;
void initScene() {
    begin1 = clock();
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent


    bones.push_back(Joint(Vector3d(0.0, 0.0, 0.0), 0.2));
    bones.push_back(Joint(bones[0].endPoint, 0.8));
    bones.push_back(Joint(bones[1].endPoint, 0.3));
    bones.push_back(Joint(bones[2].endPoint, 0.7));

    arm = System(bones);

    currGoalIndex = 0;

    myReshape(viewport.w, viewport.h);
}


clock_t now;
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50.0*zoomFactor, viewport.w/viewport.h, 1.0, 20.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    gluLookAt(  camera[0], camera[1], camera[2],
            0.0f, 1.0f,  0.0f,
            0.0f, 1.0f,  0.0f);


    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    GLfloat lightPos0[] = { 4.0f, 4.0f, 4.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

    glShadeModel(GL_FLAT);

    glPolygonMode(GL_FRONT, GL_FILL);

    now = clock();

    // double t2 = (now - begin1)%800000;
    double t2 = now - begin1;
    // cout << t2 / CLOCKS_PER_SEC << endl;

    double j = 0;
    GLfloat lightColor3[] = { 1.0f, 0.0f, 1.0f, 1.0f }; //Color (0.5, 0.5, 0.5)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor3);


    glBegin(GL_LINE_STRIP);

    while (j < 10.0){

        glVertex3d(.7 + .48*pow(sin(j), 3), 1.0 +  .39*cos(j) - 0.15 * cos(j *2) - 0.06 * cos(j *3) - 0.03 * cos(j * 4) , -1.0);
        glVertex3d(.7 + .48*pow(sin((j+.01)), 3), 1.0 +  .39*cos((j+.01)) - 0.15 * cos((j+.01) *2) - 0.06 * cos((j+.01) *3) - 0.03 * cos((j+.01) * 4) , -1.0);

        // glVertex3d(0.8 + 0.4*cos(j), 1 + 0.4*sin(2 * j), -1);
        // glVertex3d(0.8 + 0.4*cos(j + .01), 1 + 0.4*sin(2 * j + .01), -1);

        // glVertex3d(-5 + j, 0, 0);
        // glVertex3d(-5 + j+.1, 0, 0);

        j = j + .01;
    }
    glEnd();

    // infinity
    // Vector3d goal1(0.8 + 0.4*cos(.00001*t2), 1 + 0.4*sin(.00001*2*t2), -1);

    // heart
    Vector3d goal1(.7 + .48*pow(sin(.00001*t2), 3), 1.0 +  .39*cos(.00001*t2) - 0.15 * cos(.00001*t2 *2) - 0.06 * cos(.00001*t2 *3) - 0.03 * cos(.00001 * t2 * 4) , -1.0);

    // line
    // Vector3d goal1(-5+  t2*.00001, 0, 0);

    arm.update(goal1);


    arm.render();

    glPushMatrix();

    GLfloat lightColor1[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor1);

    GLfloat lightlol[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightlol);


    glBegin(GL_POLYGON);
        glVertex3d(-2.0, 0.0, -2.0);
        glVertex3d(-2.0, 0.0, 2.0);
        glVertex3d(2.0, 0.0, 2.0);
        glVertex3d(2.0, 0.0, -2.0);
    glEnd();

    glPopMatrix();

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
    glutKeyboardFunc(myKeyboard);
    glutSpecialFunc(processSpecialKeys);
    glutDisplayFunc(display);                    // function to run when its time to draw something
    glutReshapeFunc(myReshape);                  // function to run when the window gets resized
    glutIdleFunc(myFrameMove);                   // function to run when not handling any other task
    glutMainLoop();                              // infinite loop that will keep drawing and resizing and whatever else

    return 0;
}
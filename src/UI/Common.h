#pragma once

#include <QAction>
#include <QApplication>
#include <QOffscreenSurface>
#include <QOpenGLContext>
#include <QThread>
#include <QtOpenGL>

#include <Eigen/Core>

#ifdef _DEBUG
#define glDebugLog() GLError(__FILE__, __LINE__)
#else
#define glDebugLog()
#endif

// �������һ���̰߳�ȫ��OpenGL������. ע��, ������������Qt�߳���ʵ����, ��ΪOpenGL�����ı��������д���. Ȼ��, �������Ŀ������κ������߳�����Ϊ��ǰ������.
class COpenGLContextManager : public QObject
{
public:
	explicit COpenGLContextManager(int openglVersion_Magor = 2, int openglVersion_Minor = 1);

	// ͨ����OpenGL�����Ĵ��䴴�����߳��ƶ�����ǰ�߳�, ��ʹ���Ϊ��ǰ������, ��ʹOpenGL�����Ŀ���
	bool MakeCurrent();

private:
	QOffscreenSurface surface;
	QOpenGLContext context;
	QThread* parentThread;
	QThread* currentThread;
	QAction* makeCurrentAction;
};

// ��ȡOpenGL�Ĵ�����Ϣ�������������׼������(stderr)
void GLError(const char* file, int line);

Eigen::Matrix4f QMatrixToEigen(const QMatrix4x4& matrix);
QMatrix4x4 EigenToQMatrix(const Eigen::Matrix4f& matrix);











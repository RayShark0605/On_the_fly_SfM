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

// 该类管理一个线程安全的OpenGL上下文. 注意, 这个类必须在主Qt线程中实例化, 因为OpenGL上下文必须在其中创建. 然后, 该上下文可以在任何其他线程中设为当前上下文.
class COpenGLContextManager : public QObject
{
public:
	explicit COpenGLContextManager(int openglVersion_Magor = 2, int openglVersion_Minor = 1);

	// 通过将OpenGL上下文从其创建的线程移动到当前线程, 并使其成为当前上下文, 以使OpenGL上下文可用
	bool MakeCurrent();

private:
	QOffscreenSurface surface;
	QOpenGLContext context;
	QThread* parentThread;
	QThread* currentThread;
	QAction* makeCurrentAction;
};

// 获取OpenGL的错误信息并将其输出到标准错误流(stderr)
void GLError(const char* file, int line);

Eigen::Matrix4f QMatrixToEigen(const QMatrix4x4& matrix);
QMatrix4x4 EigenToQMatrix(const Eigen::Matrix4f& matrix);











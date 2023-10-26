#include "Common.h"
#include <string>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <any>
#include <Windows.h>
#include <math.h>

#include <boost/stacktrace.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

using namespace std;

// 这里之所以把Base.h中的CHECK单独复制过来而不直接include "Base.h", 是因为Base.h中include了tbb库, 而某些版本的Qt与tbb一起编译会出现编译时错误
// 参考: https://stackoverflow.com/questions/65950518/linking-intel-tbb-with-qt-compiler-msvc
#define Check(condition, ...) CheckFunc((condition), (__FILE__), (__LINE__), __VA_ARGS__)
inline void CheckFunc(bool condition, const char* file, int line, string info = "Check failed!")
{
	if (!condition)
	{
		try
		{
			auto now = chrono::system_clock::now();
			auto in_time_t = chrono::system_clock::to_time_t(now);
			stringstream ss;
			ss << put_time(localtime(&in_time_t), "%Y-%m-%d %X");

			string output = (boost::format("[%1%]%2% %3%: %4%") % ss.str() % string(file) % to_string(line) % info).str();
			throw runtime_error(output);
		}
		catch (const exception& e)
		{
			string exception = "Exception: " + string(e.what());
			OutputDebugStringA(exception.c_str());

			boost::stacktrace::stacktrace st;
			ostringstream oss;
			oss << st;
			string stacktrace = oss.str();
			OutputDebugStringA(stacktrace.c_str());

			string messageBoxContent = exception + "\n" + stacktrace;
			const size_t len = messageBoxContent.length() + 1;
			HGLOBAL hMem = GlobalAlloc(GMEM_MOVEABLE, len);
			memcpy(GlobalLock(hMem), messageBoxContent.c_str(), len);
			GlobalUnlock(hMem);
			OpenClipboard(0);
			EmptyClipboard();
			SetClipboardData(CF_TEXT, hMem);
			CloseClipboard();

			messageBoxContent = "错误信息已复制到粘贴板!\n" + messageBoxContent;
			MessageBoxA(NULL, messageBoxContent.c_str(), "Runtime Error", MB_OK | MB_ICONERROR);

			exit(EXIT_FAILURE);
		}
	}
}

COpenGLContextManager::COpenGLContextManager(int openglVersion_Magor, int openglVersion_Minor)
{
	parentThread = QThread::currentThread();
	currentThread = nullptr;
	makeCurrentAction = new QAction(this);

	Check(QCoreApplication::instance());
	Check(QCoreApplication::instance()->thread() == QThread::currentThread());

	QSurfaceFormat format;
	format.setDepthBufferSize(24);
	format.setMajorVersion(openglVersion_Magor);
	format.setMinorVersion(openglVersion_Minor);
	format.setSamples(4);
	format.setProfile(QSurfaceFormat::CompatibilityProfile);
	context.setFormat(format);

	surface.create();
	Check(context.create());
	context.makeCurrent(&surface);
	Check(context.isValid());

	connect(makeCurrentAction, &QAction::triggered, this, [this]() {
		Check(currentThread);
		context.doneCurrent();
		context.moveToThread(currentThread);
		}, Qt::BlockingQueuedConnection);
}
bool COpenGLContextManager::MakeCurrent()
{
	currentThread = QThread::currentThread();
	makeCurrentAction->trigger();
	context.makeCurrent(&surface);
	return context.isValid();
}

void GLError(const char* file, int line)
{
	GLenum errorCode(glGetError());
	while (errorCode != GL_NO_ERROR)
	{
		string errorName;
		switch (errorCode)
		{
		case GL_INVALID_OPERATION:
			errorName = "Invalid Operation";
			break;
		case GL_INVALID_ENUM:
			errorName = "Invalid Enum";
			break;
		case GL_INVALID_VALUE:
			errorName = "Invalid Value";
			break;
		case GL_OUT_OF_MEMORY:
			errorName = "Out of Memory";
			break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			errorName = "Invalid Framebuffer Operation";
			break;
		default:
			errorName = "Unknown Error";
			break;
		}
		fprintf(stderr, "OpenGL error [%s, line %i]: GL_%s", file, line, errorName.c_str());
		errorCode = glGetError();
	}
}

Eigen::Matrix4f QMatrixToEigen(const QMatrix4x4& matrix) 
{
	Eigen::Matrix4f eigenMatrix;
	for (size_t r = 0; r < 4; r++) 
	{
		for (size_t c = 0; c < 4; c++) 
		{
			eigenMatrix(r, c) = matrix(r, c);
		}
	}
	return eigenMatrix;
}
QMatrix4x4 EigenToQMatrix(const Eigen::Matrix4f& matrix)
{
	QMatrix4x4 qtMatrix;
	for (size_t r = 0; r < 4; r++)
	{
		for (size_t c = 0; c < 4; c++)
		{
			qtMatrix(r, c) = matrix(r, c);
		}
	}
	return qtMatrix;
}





#pragma once

#include <QOpenGLFunctions_3_2_Core>
#include <QtCore>
#include <QtOpenGL>

#include "../Base/Options.h"
#include "../Scene/Image.h"
#include "../Scene/Camera.h"
#include "../Scene/Database.h"
#include "Painter.h"
#include "Common.h"
#include "Colormap.h"


class CModelViewerWidget final: public QOpenGLWidget, protected QOpenGLFunctions_3_2_Core
{
public:
	CModelViewerWidget(QWidget* parent, COptions* options, CDatabase* database);
	void ShowImage(size_t imageID, size_t modelID, const Eigen::Vector4f& planeColor = Eigen::Vector4f(1, 0.1, 0, 0.6), const Eigen::Vector4f& frameColor = Eigen::Vector4f(0.8, 0.1, 0, 1));
	void ShowImages(const std::vector<size_t>& imageIDs, size_t modelID, const Eigen::Vector4f& planeColor = Eigen::Vector4f(1, 0.1, 0, 0.6), const Eigen::Vector4f& frameColor = Eigen::Vector4f(0.8, 0.1, 0, 1));




private:
	const COptions* const options;
	const CDatabase* database;
	float focusDistance;

	float point3DSize;
	float imageSize;
	float nearPlane; // 近裁剪平面
	float farPlane;
	float fieldOfView;

	float minFocusDistance;
	float maxFocusDistance;
	float focusSpeed;
	float minNearPlane;
	float maxNearPlane;
	float nearPlaneScaleSpeed;
	float minPointSize;
	float maxPointSize;
	float pointScaleSpeed;
	float minImageSize;
	float maxImageSize;
	float imageScaleSpeed;

	Eigen::Vector4f gridColor;
	Eigen::Vector4f XAxisColor; // 坐标轴的颜色
	Eigen::Vector4f YAxisColor;
	Eigen::Vector4f ZAxisColor;
	Eigen::Vector4f imagePlaneColor;
	Eigen::Vector4f imageFrameColor;

	QMatrix4x4 modelViewMatrix;
	QMatrix4x4 projectionMatrix;

	CPointPainter pointPainter;
	CLinePainter imageLinePainter;
	CTrianglePainter imageTrianglePainter;
	CLinePainter coordinateAxesPainter;
	CLinePainter coordinateGridPainter;

	CPointColormap* pointColormap;
	CImageColormap* imageColormap;


	bool isMousePressed;
	QTimer mousePressTimer;
	QPoint preMousePos;

	std::vector<std::pair<size_t, char>> selectionBuffer;
	size_t selectedImageID;
	size_t selectedPoint3DID;

	void ChangeFocusDistance(float delta); // 91
	void ChangeNearPlane(float delta);
	void ChangePointSize(float delta);
	void ChangeCameraSize(float delta);    // 94

	void RotateView(float x, float y, float preX, float preY); // 96
	void TranslateView(float x, float y, float preX, float preY); // 97
	void SelectObject(int x, int y);      // 104, TODO

	void ShowPointInfo(size_t point3DID);   // 110, TODO
	void ShowImageInfo(size_t imageID);     // 111, TODO

	void initializeGL() override;
	void resizeGL(int width, int height) override;
	void paintGL() override;                // 132
	void mousePressEvent(QMouseEvent* event) override;   
	void mouseReleaseEvent(QMouseEvent* event) override;
	void mouseMoveEvent(QMouseEvent* event) override;
	void wheelEvent(QWheelEvent* event) override;

	void Upload();                          // 143, TODO
	void UploadCoordinateGridData();        // 144, TODO
	void UploadPointData(bool isSelectionMode = false); // TODO
	void UploadPointConnectionData(); // TODO
	void UploadImageData(bool isSelectionMode = false); // TODO
	void UploadImageConnectionData();       // 148



	void ComposeProjectionMatrix();         // 151
	float GetZoomScale() const;             // 153
	float GetAspectRatio() const;           // 154
	float OrthographicWindowExtent() const; // 155

	Eigen::Vector3f PositionToArcballVector(float x, float y) const; // 157

};

















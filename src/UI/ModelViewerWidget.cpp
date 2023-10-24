#include "ModelViewerWidget.h"
#include "Common.h"

using namespace std;

#define SELECTION_BUFFER_IMAGE_IDX 0
#define SELECTION_BUFFER_POINT_IDX 1

// 从RGB颜色中生成唯一索引, 范围在[0, 256^3]之间
size_t RGBToIndex(const uint8_t r, const uint8_t g, const uint8_t b)
{
	return static_cast<size_t>(r) + static_cast<size_t>(g) * 256 + static_cast<size_t>(b) * 65536;
}
// 从由RGBToIndex生成的唯一索引中获取颜色
Eigen::Vector4f IndexToRGB(const size_t index) 
{
	Eigen::Vector4f color;
	color(0) = ((index & 0x000000FF) >> 0) / 255.0f;
	color(1) = ((index & 0x0000FF00) >> 8) / 255.0f;
	color(2) = ((index & 0x00FF0000) >> 16) / 255.0f;
	color(3) = 1.0f;
	return color;
}
// 构建影像模型用来绘制
void BuildImageModel(const CImage& image, size_t modelID, const CCamera& camera, float imageSize, const Eigen::Vector4f& planeColor, const Eigen::Vector4f& frameColor, vector<CTrianglePainterData>& triangleData, vector<CLinePainterData>& lineData)
{
	// 在OpenGL世界坐标系中生成相机维度
	const float baseCameraWidth = 1024;
	const float imageWidth = imageSize * camera.GetWidth() / baseCameraWidth;
	const float imageHeight = imageWidth * camera.GetHeight() / camera.GetWidth();
	const float imageExtentSize = max(imageWidth, imageHeight);
	const float cameraExtentSize = max(camera.GetWidth(), camera.GetHeight());
	const float cameraExtentSizeNormalized = camera.ImageToCameraThreshold(cameraExtentSize);
	const float focalLength = 2.0 * imageExtentSize / cameraExtentSizeNormalized;
	const Eigen::Matrix<float, 3, 4> invProjectionMatrix = image.GetWorldToCamera(modelID).Inverse().ToMatrix().cast<float>();

	// 投影中心, 左上角, 右上角, 右下角, 左下角
	const Eigen::Vector3f projectionCenter = invProjectionMatrix.rightCols<1>();
	const Eigen::Vector3f topLeft = invProjectionMatrix * Eigen::Vector4f(-imageWidth, imageHeight, focalLength, 1);
	const Eigen::Vector3f topRight = invProjectionMatrix * Eigen::Vector4f(imageWidth, imageHeight, focalLength, 1);
	const Eigen::Vector3f bottomRight = invProjectionMatrix * Eigen::Vector4f(imageWidth, -imageHeight, focalLength, 1);
	const Eigen::Vector3f bottomLeft = invProjectionMatrix * Eigen::Vector4f(-imageWidth, -imageHeight, focalLength, 1);


	triangleData.emplace_back(CPointPainterData(topLeft(0), topLeft(1), topLeft(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)),
		CPointPainterData(topRight(0), topRight(1), topRight(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)),
		CPointPainterData(bottomLeft(0), bottomLeft(1), bottomLeft(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)));
	triangleData.emplace_back(CPointPainterData(bottomLeft(0), bottomLeft(1), bottomLeft(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)),
		CPointPainterData(topRight(0), topRight(1), topRight(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)),
		CPointPainterData(bottomRight(0), bottomRight(1), bottomRight(2), planeColor(0), planeColor(1), planeColor(2), planeColor(3)));
	lineData.emplace_back(CPointPainterData(projectionCenter(0), projectionCenter(1), projectionCenter(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(topLeft(0), topLeft(1), topLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(projectionCenter(0), projectionCenter(1), projectionCenter(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(topRight(0), topRight(1), topRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(projectionCenter(0), projectionCenter(1), projectionCenter(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(bottomRight(0), bottomRight(1), bottomRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(projectionCenter(0), projectionCenter(1), projectionCenter(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(bottomLeft(0), bottomLeft(1), bottomLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(topLeft(0), topLeft(1), topLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(topRight(0), topRight(1), topRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(topRight(0), topRight(1), topRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(bottomRight(0), bottomRight(1), bottomRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(bottomRight(0), bottomRight(1), bottomRight(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(bottomLeft(0), bottomLeft(1), bottomLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
	lineData.emplace_back(CPointPainterData(bottomLeft(0), bottomLeft(1), bottomLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)),
		CPointPainterData(topLeft(0), topLeft(1), topLeft(2), frameColor(0), frameColor(1), frameColor(2), frameColor(3)));
}




CModelViewerWidget::CModelViewerWidget(QWidget* parent, COptions* options, CDatabase* database) : QOpenGLWidget(parent), options(options), database(database)
{
	CHECK(options && database);
	options->Check();
	isMousePressed = false;
	focusDistance = 100;
	selectedImageID = numeric_limits<size_t>::max();
	selectedPoint3DID = numeric_limits<size_t>::max();
	nearPlane = 1;
	farPlane = 100000;
	point3DSize = devicePixelRatio() * 1.0;
	imageSize = devicePixelRatio() * 0.2;
	focusDistance = 100;
	fieldOfView = 25;
	focusSpeed = 2;
	minFocusDistance = 1e-5;
	maxFocusDistance = 1e8;
	minNearPlane = 1e-3;
	maxNearPlane = 1e5;
	nearPlaneScaleSpeed = 0.02;
	minPointSize = 0.5;
	maxPointSize = 100;
	pointScaleSpeed = 0.1;
	minImageSize = 1e-6;
	maxImageSize = 1000;
	imageScaleSpeed = 0.1;

	modelViewMatrix.setToIdentity();
	modelViewMatrix.translate(0, 0, -focusDistance);
	modelViewMatrix.rotate(225, 1, 0, 0);
	modelViewMatrix.rotate(-45, 0, 1, 0);

	gridColor = Eigen::Vector4f(0.2, 0.2, 0.2, 0.6);
	XAxisColor = Eigen::Vector4f(0.9, 0, 0, 0.5);
	YAxisColor = Eigen::Vector4f(0, 0.9, 0, 0.5);
	ZAxisColor = Eigen::Vector4f(0, 0, 0.9, 0.5);
	imagePlaneColor = Eigen::Vector4f(1, 0.1, 0, 0.6);
	imageFrameColor = Eigen::Vector4f(0.8, 0.1, 0, 1);


	QSurfaceFormat format;
	format.setDepthBufferSize(24);
	format.setMajorVersion(3);
	format.setMinorVersion(2);
	format.setSamples(4);
	format.setProfile(QSurfaceFormat::CoreProfile);
#ifdef _DEBUG
	format.setOption(QSurfaceFormat::DebugContext);
#endif
	setFormat(format);
	QSurfaceFormat::setDefaultFormat(format);
}
void CModelViewerWidget::ShowImage(size_t imageID, size_t modelID, const Eigen::Vector4f& planeColor, const Eigen::Vector4f& frameColor)
{
	CHECK(database && database->IsImageExists(imageID));
	makeCurrent();
	vector<CLinePainterData> lineData;
	vector<CTrianglePainterData> tringleData;
	lineData.reserve(8);
	tringleData.reserve(2);
	
	const CImage& image = database->GetImage(imageID);
	const size_t cameraID = image.GetCameraID();
	CHECK(database->IsCameraExists(cameraID));
	const CCamera& camera = database->GetCamera(cameraID);

	BuildImageModel(image, modelID, camera, imageSize, planeColor, frameColor, tringleData, lineData);
	imageLinePainter.Upload(lineData);
	imageTrianglePainter.Upload(tringleData);
}
void CModelViewerWidget::ShowImages(const std::vector<size_t>& imageIDs, size_t modelID, const Eigen::Vector4f& planeColor, const Eigen::Vector4f& frameColor)
{
	makeCurrent();
	const size_t numImages = imageIDs.size();
	vector<CLinePainterData> lineData;
	vector<CTrianglePainterData> tringleData;
	lineData.reserve(8 * numImages);
	tringleData.reserve(2 * numImages);
	for (size_t i = 0; i < numImages; i++)
	{
		const CImage& image = database->GetImage(imageIDs[i]);
		const size_t cameraID = image.GetCameraID();
		CHECK(database->IsCameraExists(cameraID));
		const CCamera& camera = database->GetCamera(cameraID);

		BuildImageModel(image, modelID, camera, imageSize, planeColor, frameColor, tringleData, lineData);
	}
	imageLinePainter.Upload(lineData);
	imageTrianglePainter.Upload(tringleData);
}


void CModelViewerWidget::ChangeFocusDistance(float delta)
{
	if (abs(delta) < 1e-6) return;
	const float preFocusDistance = focusDistance;
	float diff = delta * GetZoomScale() * focusSpeed;
	focusDistance -= diff;
	if (focusDistance < minFocusDistance)
	{
		focusDistance = minFocusDistance;
		diff = preFocusDistance - focusDistance;
	}
	else if (focusDistance > maxFocusDistance)
	{
		focusDistance = maxFocusDistance;
		diff = preFocusDistance - focusDistance;
	}
	const Eigen::Matrix4f vmMat = QMatrixToEigen(modelViewMatrix).inverse();
	const Eigen::Vector3f tvec(0, 0, diff);
	const Eigen::Vector3f tvecRot = vmMat.block<3, 3>(0, 0) * tvec;
	modelViewMatrix.translate(tvecRot(0), tvecRot(1), tvecRot(2));
	ComposeProjectionMatrix();
	UploadCoordinateGridData();
	update();
}
void CModelViewerWidget::ChangeNearPlane(float delta)
{
	if (abs(delta) < 1e-6) return;
	nearPlane *= (1.0f + delta / 100.0f * nearPlaneScaleSpeed);
	nearPlane = max(minNearPlane, min(maxNearPlane, nearPlane));
	ComposeProjectionMatrix();
	UploadCoordinateGridData();
	update();
}
void CModelViewerWidget::ChangePointSize(float delta)
{
	if (abs(delta) < 1e-6) return;
	point3DSize *= (1.0f + delta / 100.0f * pointScaleSpeed);
	point3DSize = max(minPointSize, min(maxPointSize, point3DSize));
	update();
}
void CModelViewerWidget::ChangeCameraSize(float delta)
{
	if (abs(delta) < 1e-6) return;
	imageSize *= (1.0f + delta / 100.0f * imageScaleSpeed);
	imageSize = max(minImageSize, min(maxImageSize, imageSize));
	UploadImageData();
	update();
}
void CModelViewerWidget::RotateView(float x, float y, float preX, float preY)
{
	if (abs(x - preX) < 1e-6 && abs(y - preY) < 1e-6)return;

	// 根据Arcball方法进行旋转.
	// 参考文献: "ARCBALL: A User Interface for Specifying Three-Dimensional Orientation Using a Mouse", Ken Shoemake.

	// 确定单位球上的Arcball向量
	const Eigen::Vector3f u = PositionToArcballVector(x, y);
	const Eigen::Vector3f v = PositionToArcballVector(preX, preY);

	// 向量之间的夹角
	const float angle = 2.0f * acos(min(1.0f, u.dot(v)));

	const float kMinAngle = 1e-3f;
	if (angle > kMinAngle)
	{
		const Eigen::Matrix4f vmMat = QMatrixToEigen(modelViewMatrix).inverse();

		// 旋转坐标轴
		Eigen::Vector3f axis = vmMat.block<3, 3>(0, 0) * v.cross(u);
		axis = axis.normalized();
		// 旋转的中心是当前焦点
		const Eigen::Vector4f rot_center = vmMat * Eigen::Vector4f(0, 0, -focusDistance, 1);
		// 首先平移至旋转中心, 然后进行旋转, 最后平移回来
		modelViewMatrix.translate(rot_center(0), rot_center(1), rot_center(2));
		modelViewMatrix.rotate(RadToDeg(angle), axis(0), axis(1), axis(2));
		modelViewMatrix.translate(-rot_center(0), -rot_center(1), -rot_center(2));
		update();
	}
}
void CModelViewerWidget::TranslateView(float x, float y, float preX, float preY)
{
	if (abs(x - preX) < 1e-6 && abs(y - preY) < 1e-6)return;

	Eigen::Vector3f tvec(x - preX, preY - y, 0);
	tvec *= GetZoomScale();

	const Eigen::Matrix4f vmMat = QMatrixToEigen(modelViewMatrix).inverse();

	const Eigen::Vector3f tvec_Rot = vmMat.block<3, 3>(0, 0) * tvec;
	modelViewMatrix.translate(tvec_Rot(0), tvec_Rot(1), tvec_Rot(2));
	update();
}
void CModelViewerWidget::SelectObject(int x, int y)
{
	makeCurrent();

	// 确保抗锯齿处理不会改变对象的颜色
	glDisable(GL_MULTISAMPLE);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// 在选择模式下上传数据(每个对象一个颜色)
	UploadImageData(true);
	UploadPointData(true);

	// 以选择模式进行渲染, 使用更大的点以提高选择准确性
	const QMatrix4x4 pmv_matrix = projectionMatrix * modelViewMatrix;
	imageTrianglePainter.Render(pmv_matrix);
	pointPainter.Render(pmv_matrix, 2 * point3DSize);

	const int scaled_x = devicePixelRatio() * x;
	const int scaled_y = devicePixelRatio() * (height() - y - 1);

	QOpenGLFramebufferObjectFormat fbo_format;
	fbo_format.setSamples(0);
	QOpenGLFramebufferObject fbo(1, 1, fbo_format);

	glBindFramebuffer(GL_READ_FRAMEBUFFER, defaultFramebufferObject());
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo.handle());
	glBlitFramebuffer(scaled_x, scaled_y, scaled_x + 1, scaled_y + 1, 0, 0, 1, 1, GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);

	fbo.bind();
	array<uint8_t, 3> color;
	glReadPixels(0, 0, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, color.data());
	fbo.release();

	const size_t index = RGBToIndex(color[0], color[1], color[2]);

	if (index < selectionBuffer.size())
	{
		const char buffer_type = selectionBuffer[index].second;
		if (buffer_type == SELECTION_BUFFER_IMAGE_IDX)
		{
			selectedImageID = selectionBuffer[index].first;
			selectedPoint3DID = numeric_limits<size_t>::max();
			ShowImageInfo(selectedImageID);
		}
		else if (buffer_type == SELECTION_BUFFER_POINT_IDX)
		{
			selectedImageID = numeric_limits<size_t>::max();
			selectedPoint3DID = selectionBuffer[index].first;
			ShowPointInfo(selectionBuffer[index].first);
		}
		else
		{
			selectedImageID = numeric_limits<size_t>::max();
			selectedPoint3DID = numeric_limits<size_t>::max();
			//image_viewer_widget_->hide();
		}
	}
	else
	{
		selectedImageID = numeric_limits<size_t>::max();
		selectedPoint3DID = numeric_limits<size_t>::max();
		//image_viewer_widget_->hide();
	}

	// 重新启用, 因为在上面临时禁用了
	glEnable(GL_MULTISAMPLE);

	selectionBuffer.clear();

	UploadPointData();
	UploadImageData();
	UploadPointConnectionData();
	UploadImageConnectionData();

	update();
}
void CModelViewerWidget::ShowPointInfo(size_t point3DID)
{

}
void CModelViewerWidget::ShowImageInfo(size_t imageID)
{

}

void CModelViewerWidget::initializeGL()
{
	initializeOpenGLFunctions();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

	makeCurrent();
	pointPainter.Setup();
	imageLinePainter.Setup();
	imageTrianglePainter.Setup();
	coordinateAxesPainter.Setup();
	coordinateGridPainter.Setup();
}
void CModelViewerWidget::resizeGL(int width, int height)
{
	glViewport(0, 0, width, height);
	ComposeProjectionMatrix();
	UploadCoordinateGridData();
}
void CModelViewerWidget::paintGL()
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	const QMatrix4x4 pmvMatrix = projectionMatrix * modelViewMatrix;
	QMatrix4x4 modelViewCenterMatrix = modelViewMatrix;

	const Eigen::Vector4f rotCenter = QMatrixToEigen(modelViewMatrix).inverse() * Eigen::Vector4f(0, 0, -focusDistance, 1);
	modelViewCenterMatrix.translate(rotCenter(0), rotCenter(1), rotCenter(2));

	const QMatrix4x4 pmvcMatrix = projectionMatrix * modelViewCenterMatrix;
	coordinateAxesPainter.Render(pmvMatrix, width(), height(), 2);
	coordinateGridPainter.Render(pmvcMatrix, width(), height(), 1);

	pointPainter.Render(pmvMatrix, point3DSize);

	imageLinePainter.Render(pmvMatrix, width(), height(), 1);
	imageTrianglePainter.Render(pmvMatrix);
}
void CModelViewerWidget::mousePressEvent(QMouseEvent* event)
{
	// 选择物体(双击)
	if (mousePressTimer.isActive())
	{
		isMousePressed = false;
		mousePressTimer.stop();
		selectionBuffer.clear();
		SelectObject(event->pos().x(), event->pos().y());
	}
	else // 单击
	{
		mousePressTimer.setSingleShot(true);
		mousePressTimer.start(250);
		isMousePressed = true;
		preMousePos = event->pos();
	}
	event->accept();
}
void CModelViewerWidget::mouseReleaseEvent(QMouseEvent* event)
{
	isMousePressed = false;
	event->accept();
}
void CModelViewerWidget::mouseMoveEvent(QMouseEvent* event)
{
	if (isMousePressed)
	{
		if (event->buttons() & Qt::RightButton || (event->buttons() & Qt::LeftButton && event->modifiers() & Qt::ControlModifier))
		{
			TranslateView(event->pos().x(), event->pos().y(), preMousePos.x(), preMousePos.y());
		}
		else if (event->buttons() & Qt::LeftButton)
		{
			RotateView(event->pos().x(), event->pos().y(), preMousePos.x(), preMousePos.y());
		}
	}
	preMousePos = event->pos();
	event->accept();
}
void CModelViewerWidget::wheelEvent(QWheelEvent* event)
{
	const float delta = event->angleDelta().x() + event->angleDelta().y();
	if (event->modifiers().testFlag(Qt::ControlModifier))
	{
		ChangePointSize(delta);
	}
	else if (event->modifiers().testFlag(Qt::AltModifier))
	{
		ChangeCameraSize(delta);
	}
	else if (event->modifiers().testFlag(Qt::ShiftModifier))
	{
		ChangeNearPlane(delta);
	}
	else
	{
		ChangeFocusDistance(delta);
	}
	event->accept();
}

void CModelViewerWidget::Upload()
{
	/*CHECK(pointColormap && imageColormap);
	pointColormap->Prepare();
	imageColormap->Prepare();

	ComposeProjectionMatrix();
	UploadPointData();
	UploadImageData();
	UploadPointConnectionData();
	UploadImageConnectionData();
	update();*/
}
void CModelViewerWidget::UploadCoordinateGridData()
{
    makeCurrent();
    const float scale = GetZoomScale();

    // 视野中心网格
    vector<CLinePainterData> gridData(3);
    gridData[0].point1 = CPointPainterData(-20 * scale, 0, 0, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	gridData[0].point2 = CPointPainterData(20 * scale, 0, 0, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	gridData[1].point1 = CPointPainterData(0, -20 * scale, 0, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	gridData[1].point2 = CPointPainterData(0, 20 * scale, 0, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	gridData[2].point1 = CPointPainterData(0, 0, -20 * scale, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	gridData[2].point2 = CPointPainterData(0, 0, 20 * scale, gridColor(0), gridColor(1), gridColor(2), gridColor(3));
	coordinateGridPainter.Upload(gridData);

    // 坐标轴
    vector<CLinePainterData> axesData(3);
    axesData[0].point1 = CPointPainterData(0, 0, 0, XAxisColor(0), XAxisColor(1), XAxisColor(2), XAxisColor(3));
    axesData[0].point2 = CPointPainterData(-50 * scale, 0, 0, XAxisColor(0), XAxisColor(1), XAxisColor(2), XAxisColor(3));
	axesData[1].point1 = CPointPainterData(0, 0, 0, YAxisColor(0), YAxisColor(1), YAxisColor(2), YAxisColor(3));
    axesData[1].point2 = CPointPainterData(0, -50 * scale, 0, YAxisColor(0), YAxisColor(1), YAxisColor(2), YAxisColor(3));
    axesData[2].point1 = CPointPainterData(0, 0, 0, ZAxisColor(0), ZAxisColor(1), ZAxisColor(2), ZAxisColor(3));
    axesData[2].point2 = CPointPainterData(0, 0, -50 * scale, ZAxisColor(0), ZAxisColor(1), ZAxisColor(2), ZAxisColor(3));
	coordinateAxesPainter.Upload(axesData);
}
void CModelViewerWidget::UploadPointData(bool isSelectionMode)
{
	
}
void CModelViewerWidget::UploadPointConnectionData()
{

}
void CModelViewerWidget::UploadImageData(bool isSelectionMode)
{

}
void CModelViewerWidget::UploadImageConnectionData()
{

}



void CModelViewerWidget::ComposeProjectionMatrix()
{
	projectionMatrix.setToIdentity();
	projectionMatrix.perspective(fieldOfView, GetAspectRatio(), nearPlane, farPlane);
}
float CModelViewerWidget::GetZoomScale() const
{
	return 2.0 * tan(DegToRad(fieldOfView) / 2.0) * abs(focusDistance) / height();
}
float CModelViewerWidget::GetAspectRatio() const
{
	return width() * 1.0 / height();
}
float CModelViewerWidget::OrthographicWindowExtent() const
{
	return tan(DegToRad(fieldOfView) / 2.0) * focusDistance;
}
Eigen::Vector3f CModelViewerWidget::PositionToArcballVector(float x, float y) const
{
	Eigen::Vector3f vec(2.0f * x / width() - 1, 1 - 2.0f * y / height(), 0.0f);
	const float norm2 = vec.squaredNorm();
	if (norm2 <= 1.0)
	{
		vec.z() = std::sqrt(1.0f - norm2);
	}
	else
	{
		vec = vec.normalized();
	}
	return vec;
}






#pragma once
#include <SDL3/SDL.h>
#include <algorithm>
#include <chrono>
#include <execution>
#include <future>
#include <mutex>
#include <thread>
#include "pbrt/math/vecmath.h"
#include "pbrt/test/check.h"

namespace pbrt {
struct Color {
  static const Uint32 RED = 0xFF000000;
  static const Uint32 GREEN = 0x00FF0000;
  static const Uint32 BLUE = 0x0000FF00;
  static const Uint32 ALPHA = 0x000000FF;

  Uint8 red;
  Uint8 green;
  Uint8 blue;
  Uint8 alpha = 255;

  // Constructors
  Color(Uint8 red, Uint8 green, Uint8 blue, Uint8 alpha)
      : red(red), green(green), blue(blue), alpha(alpha)
  {}
  Color(Uint8 red, Uint8 green, Uint8 blue)
      : Color(red, green, blue, 255)
  {}
  Color(Vector3f color)
      : Color(
            std::min(255 * std::max(color.x, 0.0), 255.0),
            std::min(255 * std::max(color.y, 0.0), 255.0),
            std::min(255 * std::max(color.z, 0.0), 255.0)
        )
  {}
  Color(Uint32 color)
      : Color(
            (color & RED) >> 24,
            (color & GREEN) >> 16,
            (color & BLUE) >> 8,
            color & ALPHA
        )
  {}

  // Operators
  operator Uint32()
  {
    return (red << 24) | (green << 16) | (blue << 8) | alpha;
  }
  Color &operator*=(double scalar)
  {
    red = std::min(int(red * scalar), 255);
    green = std::min(int(green * scalar), 255);
    blue = std::min(int(blue * scalar), 255);
    alpha = std::min(int(alpha * scalar), 255);
    return *this;
  }
  Color operator*(double scalar)
  {
    Color color(red, green, blue, alpha);
    color *= scalar;
    return color;
  }
  Color &operator+=(Color &other)
  {
    red = std::min(int(red) + int(other.red), 255);
    green = std::min(int(green) + int(other.green), 255);
    blue = std::min(int(blue) + int(other.blue), 255);
    alpha = std::min(int(alpha) + int(other.alpha), 255);
    return *this;
  }
  Color operator+(Color &other)
  {
    Color color(red, green, blue, alpha);
    color += other;
    return color;
  }
  bool operator==(Color &other)
  {
    return red == other.red && green == other.green &&
           blue == other.blue && alpha == other.alpha;
  }
  friend std::ostream &operator<<(
      std::ostream &os, const Color &obj
  )
  {
    os << "[ " << int(obj.red) << ", " << int(obj.green) << ", "
       << int(obj.blue) << " ]";
    return os;
  }
};

class DisplayWindow {
 public:
  const int width;
  const int height;
  int refreshRate = 60;

  int isOpen = false;

 private:
  Uint32 *imageData;
  std::mutex windowMutex;

  bool SDL_is_initialized = false;
  SDL_Window *window;
  SDL_Renderer *renderer;
  SDL_Texture *texture;

  bool initSDL()
  {
    SDL_is_initialized = true;
    return SDL_Init(SDL_INIT_VIDEO);
  }
  bool initWindow()
  {
    window = SDL_CreateWindow(
        "Renderer", width, height, SDL_WINDOW_OPENGL
    );
    return window != nullptr;
  }
  bool initRenderer()
  {
    renderer = SDL_CreateRenderer(window, NULL);
    return renderer != nullptr;
  }
  bool initTexture()
  {
    texture = SDL_CreateTexture(
        renderer, SDL_PIXELFORMAT_RGBA8888,
        SDL_TEXTUREACCESS_STREAMING, width, height
    );
    return texture != nullptr;
  }
  bool initAll()
  {
    if (!(ASSERT(
              initSDL(),
              "Error Initializing SDL: " << SDL_GetError()
          ) &&
          ASSERT(
              initWindow(),
              "Error Creating window: " << SDL_GetError()
          ) &&
          ASSERT(
              initRenderer(),
              "Error Creating renderer: " << SDL_GetError()
          ) &&
          ASSERT(
              initTexture(),
              "Error creating texture: " << SDL_GetError()
          )))
    {
      ERROR("Failed to initialize DisplayWindow.  Destroying.");
      DisplayWindow::~DisplayWindow();
      return false;
    }
    else {
      return true;
    }
  }

  bool updateTexture()
  {
    if (!ASSERT(
            SDL_UpdateTexture(
                texture, NULL, imageData, width * sizeof(Uint32)
            ),
            "Failed to update texture: " << SDL_GetError()
        ))
    {
      DisplayWindow::~DisplayWindow();
      return false;
    }
    return true;
  }

  std::promise<void>
      windowClosed;  // Signals when window closes
  std::future<void> waitClosed;
  std::mutex closingWindow;
  void closeWindow()
  {
    closingWindow.lock();

    if (!isOpen) {
      closingWindow.unlock();
      return;
    }

    isOpen = false;
    windowClosed.set_value();
    if (window != nullptr) SDL_DestroyWindow(window);

    closingWindow.unlock();
  }

  void Draw()
  {
    // Copy the texture to the renderer
    SDL_RenderTextureRotated(
        renderer, texture, NULL, NULL, 0, NULL,
        SDL_FLIP_VERTICAL
    );

    // Present the frame
    SDL_RenderPresent(renderer);
  }

 public:
  DisplayWindow(
      const int width, const int height, int refreshRate
  )
      : width(width), height(height), refreshRate(refreshRate)
  {
    imageData = new Uint32[width * height] {0x000000FF};
  }
  DisplayWindow(const int width, const int height)
      : DisplayWindow(width, height, 60)
  {}

  ~DisplayWindow()
  {
    closeWindow();
    if (renderer != nullptr) SDL_DestroyRenderer(renderer);
    if (texture != nullptr) SDL_DestroyTexture(texture);
    if (SDL_is_initialized) SDL_Quit();
    if (imageData != nullptr) delete[] imageData;
  }
  void Open()
  {
    windowClosed = std::promise<void>();
    waitClosed = windowClosed.get_future();
    isOpen = true;

    std::thread windowThread([this]() {
      std::lock_guard<std::mutex> lock(windowMutex);
      initAll();
      while (true) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
          if (event.type == SDL_EVENT_QUIT) {
            closeWindow();
            return;
          }
        }

        updateTexture();
        Draw();

        // Wait before refreshing
        std::this_thread::sleep_for(
            std::chrono::milliseconds(int(1000.0 / refreshRate))
        );
      }
    });

    windowThread.detach();
  }
  void WaitUntilClosed() { waitClosed.get(); }
  bool UpdateImage(
      const Uint32 *newImage, const unsigned int size
  )
  {
    if (size > width * height) return false;

    std::copy_n(std::execution::par, newImage, size, imageData);
    return updateTexture();
  }
  // newImage is an array of Vector3f in range [0,1]
  // Values larger than 1 are clamped to 1
  bool UpdateImage(
      const Vector3f *newImage, const unsigned int size
  )
  {
    if (size > width * height) return false;

    std::transform(
        std::execution::par, newImage, newImage + size,
        imageData, [](const Vector3f &v) { return Color(v); }
    );

    return updateTexture();
  }
  void SetPixel(const unsigned int index, const Uint32 color)
  {
    if (index >= width * height) return;
    imageData[index] = color;
  }
};
}  // namespace pbrt

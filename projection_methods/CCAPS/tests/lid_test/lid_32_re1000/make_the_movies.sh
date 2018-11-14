ffmpeg -framerate 10 -i "stream_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" stream.mov

ffmpeg -framerate 10 -i "ghia_u_cent_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" ghia.mov

ffmpeg -framerate 10 -i "u_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" u.mov

ffmpeg -framerate 10 -i "v_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" v.mov

ffmpeg -framerate 10 -i "div_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" div.mov

ffmpeg -framerate 10 -i "curl_%04d.png" -c:v libx264 -crf 23  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" curl.mov




